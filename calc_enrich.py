import warnings
import MDAnalysis as mda
import pandas as pd
import numpy as np
import sys

import lipyphilic
from lipyphilic.lib.neighbours import Neighbours



def main(u,lipids,frame0):

    # filter warnings
    warnings.filterwarnings("ignore")


    # BEAD SELECTION (too long if use resname)
    sel_str = "name NP ROH"


    # calculate whether bead is within the cutoff for listed bead names
    # results available in neighbours.Neighbours as np array of matrices
    print('\nCalculating neighbours:')
    neighbours = Neighbours(universe=u,
                            lipid_sel=sel_str,
                            cutoff=12.0
                            )
    neighbours.run(start=frame0,stop=None,step=None,verbose=True)

    # calculate enrichment index of lipid species
    # the number of each neighbour species B around a given reference species A
    # is normalized by the average number of species B around any lipid.
    counts, enrichment = neighbours.count_neighbours(return_enrichment=True)

    # enrichment is a n>2 column df where a given component is calculated against another
    # e.g. Label = 'CHOL' has columns 'feCHOL' and 'feCIL'
    # where it's enrichment/depletion of chol with respect to chol or CIL
    # the frames column gives the frame of the calculation
    # this is for beads (head) within the cutoff range (12 A default)
    # print(enrichment)

    # using labels with fe as prefix, this stands for fractional enrichment
    # each df will have columns: Label, frame, fe..., fe..., ...
    for lipid in lipids:
        df = enrichment.loc[enrichment['Label']==lipid]
        print(f'\nLipid: {lipid}')

        # within the dataframe iterate through the columns to get all values
        for LIPID in lipids:

            # write column header
            column = 'fe' + LIPID

            # access column from df
            series = df.loc[:,column]

            # calculate the mean fractional enrichment/depletion
            mean = float(series.describe().loc[['mean']])

            # calculate the standard error of the dataframe
            se = float(series.sem())

            # print statistics to terminal
            print(f'Average fe with respect to {LIPID} = {mean:.2f}')

            # print("\n\n%s" %series.describe())

    return enrichment
