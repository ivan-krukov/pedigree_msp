import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import argparse


class PedFiller:
    def __init__(self, ped_df):
        self.ped_df = ped_df
        self.min_ID = 0
        self.probands = set()
        self.new_ped = None
        # TODO: Would be good to have 'parent_cohorts' where
        # we can choose from pre-existing couples rather than
        # just random individuals from each cohort
        # - could build cohorts by collecting parents who had a
        #   child at the given time, that way the generation time
        #   is chosen automatically from the pedigree itself.
        self.node_cohorts = defaultdict(list)
        self.founder_cohorts = defaultdict(list)

    @staticmethod
    def read_ped(pedfile, header=0, sep=' ', dtype=int, usecols=(0, 1, 2)):
        df = pd.read_csv(pedfile,
                         header=header,
                         sep=sep,
                         dtype=dtype,
                         usecols=usecols)
        for expected_col in ['ind', 'mother', 'father']:
            assert(expected_col in df.columns)

        return PedFiller(df)

    def partial_order_times(self):
        ind_dict = dict(zip(self.ped_df.ind, range(len(self.ped_df.ind))))
        not_mothers = set(self.ped_df['ind']).difference(self.ped_df['mother'])
        self.probands = not_mothers.difference(self.ped_df['father'])

        climbers = self.probands
        times = np.zeros(len(self.ped_df['ind']))
        t = 0
        while len(climbers) > 0:
            next_climbers = set()
            for c in tqdm(climbers,
                          leave=False,
                          desc="Setting time {}".format(t)):
                idx = ind_dict[c]
                time = times[idx]
                if t > time:
                    times[idx] = t

                mother = self.ped_df['mother'][idx]
                father = self.ped_df['father'][idx]
                if mother != 0:
                    next_climbers.add(mother)
                if father != 0:
                    next_climbers.add(father)

            climbers = set(next_climbers)
            t += 1

        self.ped_df['time'] = times

    def build_cohorts(self):
        """
        Collects individuals into cohorts of identical times
        """
        for i, row in tqdm(self.ped_df.iterrows(),
                           total=self.ped_df.shape[0],
                           desc='Building cohorts'):
            if row.mother == 0 and row.father == 0:
                self.founder_cohorts[row.time].append(int(row.ind))

            # Node cohorts contain both founders and internal nodes
            self.node_cohorts[row.time].append((row.ind))

    def get_wf_parents(self, n, N, min_ID, monogamous=True):
        """
        Draws new parents for unconnected recent founders from a
        randomly-mating Wright-Fisher population
        """
        if monogamous is not True:
            raise NotImplementedError

        num_couples = N / 2
        mothers = np.arange(min_ID,
                            min_ID + num_couples,
                            dtype=int
                            ).reshape(-1, 1)
        fathers = np.arange(min_ID + num_couples,
                            min_ID + (num_couples * 2),
                            dtype=int
                            ).reshape(-1, 1)

        if monogamous is True:
            np.random.shuffle(fathers)
            couple_indices = np.random.randint(0, num_couples, size=n)
            couples = np.concatenate([mothers, fathers],
                                     axis=1)[couple_indices]

        return couples

    def complete_pedigree(self, N=100, reconnection_rate=0.05, max_gen_diff=1):
        """
        Draws new connections for recent founders from parents within
        max_gen_diff of the recent founder's time.
        """
        assert(len(self.founder_cohorts) > 0)
        assert(len(self.node_cohorts) > 0)

        max_time = int(np.max(self.ped_df['time']))
        max_ID = np.max(self.ped_df['ind'])

        self.new_ped_rows = []
        for time in tqdm(range(max_time), desc='Drawing new connections'):
            founders = self.founder_cohorts[time]
            if len(set(founders)) != len(founders):
                print("Error at time", time)
                print(founders)
                break
            if len(founders) == 0:
                continue

            possible_parents = []
            for i in range(1, max_gen_diff + 1):
                possible_parents.extend(self.node_cohorts[time + i])
            assert(len(possible_parents) == len(set(possible_parents)))

            num_to_reconnect = np.random.binomial(len(founders),
                                                  reconnection_rate)
            np.random.shuffle(founders)
            to_reconnect = founders[:num_to_reconnect]
            not_reconnected = founders[num_to_reconnect:]

            new_parents = self.get_wf_parents(len(not_reconnected), N, max_ID)
            max_ID = np.max(new_parents)
            for founder, parents in zip(not_reconnected, new_parents):
                mother, father = parents
                row = [founder, mother, father, time]
                self.new_ped_rows.append(row)

            # Add new parents to next founder generation
            self.founder_cohorts[time + 1].extend(
                set(list(new_parents.ravel()))
            )

            for founder in to_reconnect:
                mother, father = np.random.choice(possible_parents, size=2)
                row = [founder, mother, father, time]
                self.new_ped_rows.append(row)

    def collect_rows(self):
        """
        Sorts original and new individuals into unchanged/updated/new
        categories, then combines them into a new array
        """
        new_df = pd.DataFrame(self.new_ped_rows,
                              columns=['ind', 'father', 'mother', 'time']
                              ).set_index('ind').sort_index()
        original_df = self.ped_df.set_index('ind').sort_index()

        self.merged_rows = []
        self.unchanged_rows = []
        self.new_rows = []
        new_iter = new_df.iterrows()
        new_founder, new_row = next(new_iter)

        # Replace reconnected founder rows
        for ind, row in tqdm(original_df.iterrows(),
                             total=original_df.shape[0],
                             desc="Merging new connections"):
            if new_founder == ind:
                self.merged_rows.append(
                    [new_founder,
                     int(new_row.father),
                     int(new_row.mother),
                     new_row.time]
                )
                new_founder, new_row = next(new_iter)
            else:
                self.unchanged_rows.append(
                    [ind,
                     int(row.father),
                     int(row.mother),
                     row.time]
                )

        # Add new individuals
        for founder, row in new_iter:
            self.new_rows.append(
                [founder,
                 int(row.father),
                 int(row.mother),
                 row.time]
            )

        assert(len(self.unchanged_rows) + len(self.merged_rows) ==
               self.ped_df.shape[0])

        self.new_ped = self.unchanged_rows + self.merged_rows + self.new_rows

    def write_filled_ped(self, path):
        with open(os.path.expanduser(path), 'w') as f:
            f.write('ind\tfather\tmother\ttime\n')
            for row in self.new_ped:
                f.write('\t'.join([str(x) for x in row]) + '\n')


def main(args):
    sep = args.separator
    if args.separator.lower() == 'tab':
        sep = '\t'
    elif args.separator.lower() == 'space':
        sep = ' '

    P = PedFiller.read_ped(args.pedfile, sep=sep)

    P.partial_order_times()
    P.build_cohorts()
    P.complete_pedigree(args.population_size, args.reconnection_rate)
    P.collect_rows()
    P.write_filled_ped(args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--pedfile', metavar='', required=True,
                        help="""Path to pedigree text file""")
    parser.add_argument('-o', '--outfile', metavar='', required=True,
                        help="""Path to save filled pedigree""")
    parser.add_argument('-s', '--separator', metavar='', default='tab',
                        help="""Separator used in 'pedfile. For <space> use
                        'space'. Default: <tab>""")
    parser.add_argument('-r', '--reconnection-rate', metavar='', default=0.05,
                        help="""Rate at which recent founders will reconnect to
                        the pedigree. Default: %(default).3f""")
    parser.add_argument('-N', '--population-size', metavar='', default=100,
                        help="""Size of randomly-mating population containing
                        founder lineages which have not reconnected. Default:
                        %(default)d""")

    args = parser.parse_args()
    main(args)
