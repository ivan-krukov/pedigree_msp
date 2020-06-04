import sys
import os
import argparse
import subprocess
import numpy as np
import tempfile

sys.path.insert(0, 'msprime')
import msprime # NOQA


class Runner:
    def __init__(self, args):
        self.test = args.test
        self.verbose = (args.verbose or args.test)

    def run(self, cmd):
        prefix = "Running:"
        if self.test:
            prefix = "Testing:"

        if self.verbose:
            print(prefix, cmd)

        if not self.test:
            subprocess.run(cmd, shell=True, check=True)

def test_afs_calculation(ts_file='cached/balsac_140.tskit'):
    if not os.path.exists(ts_file):
        
        ped = msprime.Pedigree.read_txt('data/balsac_140.tsv')
        ped.set_samples(num_samples=14)
        des = [msprime.SimulationModelChange(max(ped.time))]
        sim = msprime.simulate(14,
                               Ne=1000,
                               pedigree=ped,
                               model='wf_ped',
                               length=1e8,
                               mutation_rate=1e-8,
                               recombination_rate=1e-8,
                               demographic_events=des)

        sim.dump('cached/balsac_140.tskit')
    
    ts = msprime.load(ts_file)
    sample_nodes = ts.samples()
    afs = np.zeros(len(sample_nodes)+1)
    for tree in ts.trees():
        for mutation in tree.mutations():
            count = tree.num_samples(mutation.node)
            frequency = count / len(sample_nodes)
            afs[count] += 1

    ts_afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
    
    assert np.allclose(afs, ts_afs)


# TODO - provide a way to subsample based on a file with IDs
def ts_to_bcf_single(ts_file, out_file, runner, af_cutoff=0, n_subsample=None, remove_singletons=False):

    ts = msprime.load(ts_file)
    if ts.num_individuals == 0:
        bcf_cmd = "tskit vcf --ploidy 2 {} | bcftools view -O b > {}".format(
                           ts_file, out_file)
        runner.run(bcf_cmd)
    else:
        # remove sites based on allele frequency cutoff
        sample_nodes = ts.samples()

        if (af_cutoff != 0) or remove_singletons:
            sites_to_delete = []
            for tree in ts.trees():
                for mutation in tree.mutations():
                    # should this be count-1 since the node itself is included?
                    count = tree.num_samples(mutation.node)
                    # note that this is non-polarized frequency ( > 50% possible)
                    if (count / len(sample_nodes)) < af_cutoff:
                        sites_to_delete.append(mutation.site)
                    if remove_singletons and count == 1:
                        sites_to_delete.append(mutation.site)

            # bulk-remove sites
            tables = ts.dump_tables()
            tables.delete_sites(sites_to_delete, record_provenance=False)
            ts = tables.tree_sequence()
            
        sample_individuals = []
        for ind in ts.individuals():
            if len(ind.nodes) == 0:
                continue
            # diploid - two nodes, both the same individual
            if ind.nodes[0] in sample_nodes:
                assert len(ind.nodes) == 2
                assert ind.nodes[1] in sample_nodes
                sample_individuals.append(ind)

        if n_subsample is not None:
            sample_individuals = np.random.choice(sample_individuals, n_subsample, replace=False)
                
        # TODO: This will probably fail if metadata isn't present
        sample_names = [ind.metadata.decode('utf8') for ind in sample_individuals]
        sample_ids = [ind.id for ind in  sample_individuals]
        
        
        read_fd, write_fd = os.pipe()
        write_pipe = os.fdopen(write_fd, "w")
        with open(out_file, "w") as f:
            proc = subprocess.Popen(
                ["bcftools", "view", "-O", "b"], stdin=read_fd, stdout=f
            )
            ts.write_vcf(
                write_pipe, individuals=sample_ids, individual_names=sample_names
            )
            write_pipe.close()
            os.close(read_fd)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(
                    "bcftools failed with status:", proc.returncode
                )


def bcf_convert_chrom(bcf_file, chrom_num, runner):
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write('1\t' + str(chrom_num) + '\n')

        convert_chrom_cmd = (
                'bcftools annotate --rename-chrs {} {} | '
                'bcftools view -O b > tmp && mv tmp {} && '
                'rm {}').format(
                f.name, bcf_file, bcf_file, f.name)
        runner.run(convert_chrom_cmd)


def concat_bcf(bcf_files, out_file, runner):
    concat_chrom_cmd = 'bcftools concat -o {} -O b --threads 2'

    for f in bcf_files:
        concat_chrom_cmd += ' ' + f

    concat_chrom_cmd = concat_chrom_cmd.format(out_file)

    runner.run(concat_chrom_cmd)


def main(args):
    out_file = os.path.expanduser(args.out_file)
    runner = Runner(args)

    tmp_bcf_files = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        for i, tsf in enumerate(args.ts_file):
            chrom_num = i + 1
            tmp_bcf_file = os.path.join(tmpdirname, '.tmp' + str(i) + '.bcf')
            tmp_bcf_files.append(tmp_bcf_file)

            ts_to_bcf_single(tsf, tmp_bcf_file, runner, af_cutoff=args.af_cutoff, n_subsample=args.n_subsample, remove_singletons=args.remove_singletons)
            bcf_convert_chrom(tmp_bcf_file, chrom_num, runner)

        concat_bcf(tmp_bcf_files, out_file, runner)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--ts-file', nargs="*", required=True)
    parser.add_argument('-o', '--out-file', required=True)
    parser.add_argument('-s', '--n_subsample', type=int, default=None, help="Randomly sumbsample this many individuals. Default - all")
    parser.add_argument('-f', '--af_cutoff', type=float, default=0, help="Drop sites below this allele frequency cutoff. Default - 0")
    parser.add_argument('-r', '--remove-singletons', action='store_true', help="Drop sites with a single mutation")
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-T', '--test', action='store_true')

    args = parser.parse_args()

    main(args)
