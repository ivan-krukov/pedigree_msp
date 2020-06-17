import msprime
import sys
from argparse import ArgumentParser

parser = ArgumentParser('sim_continue')

parser.add_argument('-i', '--input-ts', help='input incomplete tree sequence')
parser.add_argument('-o', '--output-ts', help='output tree sequence file')
parser.add_argument('-N', '--population-size', help='population size', type=int, default=100)

args = parser.parse_args()

ts = msprime.load(args.input_ts)

sim = msprime.simulate(from_ts = ts, model='dtwf', Ne=args.population_size)

sim.dump(args.output_ts)
