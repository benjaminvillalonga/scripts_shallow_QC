"""Scrpit to generate random circuits according to Sergio Boixo's prescription.

Parameters:
  --L1: first lateral dimension of the circuit.
  --L2: second lateral dimension of the circuit.
  --depth: depth of the circuit (greater than 1).
  --seed: random seed to generate the circuit.
  --outfile: output file for the circuit.

Raises:
  ValueError: if the depth is smaller than 2.

Example:
  ``python generate_circuit.py --L1 7 --L2 7 --depth 100 --seed 0 --outfile example_circuits/inst_7_7_100_0``

  This command will generate a file *inst_7_7_100_0* in folder *example_circuits/* of a circuits with dimensions 7x7 and depth 100, randomly generated from seed 0.

Creates a file *outfile* (no extension added) with the circuit in the format:
  * First line has the number of qubits.
  * Subsequent lines have either 3 or 4 fields: cycle, gate (*h*, *x_1_2*, *y_1_2*, *t*, or *cz*).

First lay down the Hadamard gates on the first cycle.  Then lay down the layers of cz gates according to the prescription.  Finnaly, fill in the *t*, *x_1_2*, and *y_1_2* gates in between.  Last step would be to print all layers in order according to the cycles.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import itertools as it
import argparse
import logging

parser = argparse.ArgumentParser(description='Generate a random circuit.')
parser.add_argument('--L1', type=int,
                    help='First lateral dimension of the circuit.')
parser.add_argument('--L2', type=int,
                    help='Second lateral dimension of the circuit.')
parser.add_argument('--depth', type=int,
                    help='Depth of the circuit (greater than 1)')
parser.add_argument('--seed', type=int,
                    help='Random seed to generate the circuit.')
parser.add_argument('--outfile', type=str,
                    help='Output filename.')

_GATES = ['h', 't', 'x_1_2', 'y_1_2', 'cz']

class Circuit(object):
  """A CircuitGrid provides the skeleton of a circuit while it is being built.

  Args:
    L1: first lateral dimension of the circuit.
    L2: second lateral dimension of the circuit.
    depth: depth of the circuit.

  Attributes:
    L1: first lateral dimension of the circuit.
    L2: second lateral dimension of the circuit.
    depth: depth of the circuit.
    grid: a str numpy.array with the sites (i1, i2, cycle) taken by a gate.
    gates_at_cycle: list of gates (cycle, gate, idx(s)) at each cycle.
  """
  def __init__(self, L1, L2, depth):
    self.L1 = L1
    self.L2 = L2
    self.depth = depth
    self.grid = np.zeros((L1, L2, depth), dtype=object)
    self.gates_at_cycle = [[] for cycle in range(depth)]
    logging.info('Created a Circuit {}'.format(self.__repr__()))

  def __repr__(self):
    """String representation of the circuit."""
    repr_cg = '<Circuit: L1={}, L2={}, depth={}>'.format(self.L1, self.L2,\
                        self.depth)
    return repr_cg

  def coords_to_idx(self, i1, i2):
    """Function to obtain the idx of a qubit from its 2D coordinates.

    Args:
      i1: first spatial coordinate.
      i2: second spatial coordinate.

    Returns:
      idx of the qubit corresponting to the coordinates given.
    """
    idx = i1 + i2*self.L1
    return idx

  def idx_to_coords(self, idx):
    """Function to obtain the 2D coordinates from the idx of a qubit.

    Args:
      idx: int index of a qubit.

    Returns:
      (i1, i2) tuple of ints with the 2D coordinates.
    """
    i1 = idx%self.L1
    i2 = idx//self.L1
    return (i1, i2)

  def set_gate(self, gate, coords, cycle):
    """Sets a gate at the position specified by *coords* and *cycle*.

    Args:
      gate: one of 'h', 't', 'x_1_2', 'y_1_2', 'cz'.
      coords: coordinates of the gate in the format (i1, i2) if it is a one-qubit gate, or (i1, i2, j1, j2) if it is a two-qubit gate.
      cycle: cycle at which the gate is applied time coordinate).

    Raises:
      ValueError: if gate doesn't belong to the valid set or if the gate and the coordinates don't match.
      ValueError: if the sites in the grid are taken by another gate.

    Todo: Check that the 'cz' gates are set on neighbors.
    """
    if (gate not in _GATES) or (gate=='cz' and len(coords)!=4) or\
      (gate!='cz' and len(coords)!=2):
      msg = 'gate has to be valid and coords match the chosen gate.'
      raise ValueError(msg)
    if len(coords)==2:
      (i1, i2) = coords
      if self.grid[i1, i2, cycle]:
        msg = 'The site where you attempt to set the gate is already taken.'
        raise ValueError(msg)
    if len(coords)==4:
      (i1, i2, j1, j2) = coords
      if self.grid[i1, i2, cycle] or self.grid[j1, j2, cycle]:
        msg = 'The sites where you attempt to set the gate are already taken.'
        raise ValueError(msg)

    logging.debug('Setting {} at coordinates {} and cycle {}'\
            .format(gate, coords, cycle))

    if len(coords)==2:
      idx = self.coords_to_idx(i1, i2)
      self.grid[i1, i2, cycle] = gate
      self.gates_at_cycle[cycle].append('{} {} {}'.format(cycle, gate, idx))

    if len(coords)==4:
      idx1 = self.coords_to_idx(i1, i2)
      idx2 = self.coords_to_idx(j1, j2)
      self.grid[i1, i2, cycle] = gate
      self.grid[j1, j2, cycle] = gate
      self.gates_at_cycle[cycle].append('{} {} {} {}'.format(cycle,\
                                          gate, idx1, idx2))

  def circuit_to_file(self, filename):
    """Print circuit to a file *filename*

    Args:
      filename: name of the output file.
    """
    file_out = open(filename, 'w')
    file_out.write('{}\n'.format(self.L1*self.L2))
    for cycle in range(self.depth):
      for line in self.gates_at_cycle[cycle]:
        file_out.write('{}\n'.format(line))


def cz_layer(L1, L2, pattern_num):
  """Function that generates the cz layer corresponding to the pattern_num.

  Args:
    L1: first spatial dimension.
    L2: second spatial dimension.
    pattern_num: number that labels the 8 different patterns of cz layers.

  Returns:
    A list with four-tuples with the (i1, i2, j1, j2) coords of the cz gates.

  Raises:
    ValueError: if pattern_num is not in [1,8].
  """
  if pattern_num<1 or pattern_num>8:
    msg = 'pattern_num must be in the range [1,8].'
    raise ValueError(msg)

  coords_list = []
  if pattern_num==1:
    for i2 in range(L2):
      if i2%2==0:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(0, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(2, L1-1, 4)]
  elif pattern_num==2:
    for i1 in range(L2):
      if i1%2==0:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(1, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(3, L1-1, 4)]
  elif pattern_num==3:
    for i2 in range(L2):
      if i2%2==0:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(1, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(3, L1-1, 4)]
  elif pattern_num==4:
    for i1 in range(L2):
      if i1%2==0:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(0, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(2, L1-1, 4)]
  elif pattern_num==5:
    for i2 in range(L2):
      if i2%2==0:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(2, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(0, L1-1, 4)]
  elif pattern_num==6:
    for i1 in range(L2):
      if i1%2==0:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(3, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(1, L1-1, 4)]
  elif pattern_num==7:
    for i2 in range(L2):
      if i2%2==0:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(1, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1+1, i2) for i1 in range(3, L1-1, 4)]
  elif pattern_num==8:
    for i1 in range(L2):
      if i1%2==0:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(2, L1-1, 4)]
      else:
        coords_list += [(i1, i2, i1, i2+1) for i2 in range(0, L1-1, 4)]
  
  return coords_list


def main(FLAGS):
  np.random.seed(FLAGS.seed)

  L1, L2, depth = FLAGS.L1, FLAGS.L2, FLAGS.depth
  num_qubits = L1*L2
  my_circuit = Circuit(L1, L2, depth)
  
  # Set first layer of h gates.
  for idx in range(num_qubits):
    (i1, i2) = my_circuit.idx_to_coords(idx)
    my_circuit.set_gate('h', (i1, i2), 0)

  # Set cz layers.
  for cycle in range(1, depth):
    layer = cz_layer(L1, L2, (cycle-1)%8+1)
    for coords in layer:
      my_circuit.set_gate('cz', coords, cycle)

  # Set one-qubit gates in between the cz ones.
  fill_in_gates = ['t', 'x_1_2', 'y_1_2']
  for idx in range(num_qubits):
    (i1, i2) = my_circuit.idx_to_coords(idx)
    coords = (i1, i2)
    world_line = my_circuit.grid[i1, i2, :]
    first_t_set = False 
    last_gate = ''
    for cycle, site_content in enumerate(world_line):
      if not first_t_set:
        if not site_content:
          my_circuit.set_gate('t', coords, cycle)
          first_t_set = True
          last_gate = 't'
      else:
        if world_line[cycle-1]=='cz' and not site_content:
          available_gates = list(fill_in_gates)
          available_gates.remove(last_gate)
          gate = available_gates[np.random.randint(len(available_gates))]
          my_circuit.set_gate(gate, coords, cycle)
          last_gate = gate
  
  # Print circuit to file.
  my_circuit.circuit_to_file(FLAGS.outfile)


if __name__ == '__main__':
  FLAGS = parser.parse_args()
  main(FLAGS)
