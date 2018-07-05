"""Script to generate a cnf file (to be read by quickbb) from a circuit.

Parameters:
  --circuit_file: name of the input (circuit) file.
  --outfile: name of the file where the graph is written.
  --max_cycle: maximum cycle included in the circuit.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import logging

import circuit_graph as circuit

parser = argparse.ArgumentParser(description='Generate graph from circuit.')
parser.add_argument('--circuit_file', type=str,
                    help='File that contains the circuit information.')
parser.add_argument('--outfile', type=str,
                    help='Output file for the graph of the circuit.')
parser.add_argument('--max_cycle', type=int,
                   help = 'Include gates up to cycle <= max_cycle.')


def main(FLAGS):

  # Read circuit.
  with open(FLAGS.circuit_file, 'r') as input_file:
    lines = input_file.readlines()

  # Create graph in Sergio Boixo's format.
  tuples, n_qubits = circuit.read_circuit_lines(iter(lines), FLAGS.max_cycle)
  circuit_graph = circuit.get_circuit_graph(tuples, n_qubits)

  # Count number of variables and cliques.
  variables = set()
  gate_strings = []
  for gate in circuit_graph.graph:
    if len(gate)==2:
      variables |= set([gate[1]])
      gate_strings.append('{} 0'.format(gate[1]))
    if len(gate)==3:
      variables |= set([gate[1], gate[2]])
      gate_strings.append('{} {} 0'.format(gate[1], gate[2]))
  num_variables = len(variables)

  # Write to file.
  file_out = open(FLAGS.outfile, 'w')
  file_out.write('p cnf {} {}\n'.format(num_variables, len(gate_strings)))
  for line in gate_strings:
    file_out.write('{}\n'.format(line))
  file_out.close()


if __name__ == '__main__':
  FLAGS = parser.parse_args()
  main(FLAGS)
