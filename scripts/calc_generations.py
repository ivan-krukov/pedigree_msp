import msprime
import sys

ts = msprime.load(sys.argv[1])
gen = max(n.time for n in ts.nodes())
print(gen)
