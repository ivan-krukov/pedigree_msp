import msprime
import sys
from argparse import ArgumentParser
import json

parser = ArgumentParser('sim_continue')

parser.add_argument('-i', '--input-ts', help='input incomplete tree sequence')
parser.add_argument('-o', '--output-ts', help='output tree sequence file')
parser.add_argument('-N', '--population-size', help='population size', type=int, default=100)

args = parser.parse_args()

ts = msprime.load(args.input_ts)

# read the recombination map from the provenance
provenance = json.loads(ts.provenance(0).record)
ts_recomb_map = provenance["parameters"]["recombination_map"]
rm = msprime.RecombinationMap(ts_recomb_map["positions"], ts_recomb_map["rates"], discrete=ts_recomb_map["discrete"])

sim = msprime.simulate(from_ts = ts, model='dtwf', Ne=args.population_size, recombination_map=rm)

sim.dump(args.output_ts)
