import msprime as msp

ped = msp.Pedigree.read_txt('data/balsac_140.tsv')
ped.set_samples(num_samples=14)
des = [msp.SimulationModelChange(max(ped.time))]

sim = msp.simulate(14,
                   Ne=1000,
                   pedigree=ped,
                   model='wf_ped',
                   length=1e8,
                   mutation_rate=1e-8,
                   recombination_rate=1e-8,
                   demographic_events=des)

sim.dump('cached/balsac_140.tskit')
