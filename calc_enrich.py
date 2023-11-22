import warnings
import MDAnalysis as mda
import pandas as pd
import numpy as np
import sys

import lipyphilic
from lipyphilic.lib.neighbours import Neighbours


def main(u,sysComp,frame0):

    # NB: only vlaid for systems with cholesterol

    # filter warnings
    warnings.filterwarnings("ignore")

    # calculate whether bead is within the cutoff for listed bead names
    neighbours = Neighbours(
        universe=u,
        lipid_sel="name N1 ROH",
        cutoff=12.0
        )

    # results available in neighbours.Neighbours as np array of matrices
    neighbours.run(start=frame0,stop=None,step=None,verbose=True)

    # calculate enrichment index of lipid species
    # the number of each neighbour species B around a given reference species A
    # is normalized by the average number of species B around any lipid.
    # returns two pandas dfs
    counts, enrichment = neighbours.count_neighbours(return_enrichment=True)

    #
    enrichment_filtered = enrichment.loc[enrichment['Label']=='CHOL']
    # print(enrichment_filtered)

    y = [enrichment_filtered.loc[:,'feCHOL'],enrichment_filtered.loc[:,'feMC3H']]
    # print(y)

    print(f'\nAverage chol enrichment = {np.mean(y[0][-4000:]):.2f}')
    print(f'Std chol enrichment = {np.std(y[0][-4000:]):.2f}')
    print(f'Average CIL enrichment = {np.mean(y[1][-4000:]):.2f}')
    print(f'Std CIL enrichment = {np.std(y[1][-4000:]):.2f}')

    # define x axis (number of frames) based on length of y points
    frames = [i for i in range(len(y[0]))]
    x = [frames,frames]

    # enrichment_chol = enrichment.loc[enrichment['Label']=='CHL1']
    # enrichment_chol.plot(title='chol',x='Frame',y=['feCHL1','feDLMC3'])#,'feADE'])
    # # print(enrichment_chol)

    # enrichment_mc3 = enrichment.loc[enrichment['Label']=='DLMC3']
    # enrichment_mc3.plot(title='mc3',x='Frame',y=['feCHL1','feDLMC3'])#,'feADE'])
    # print(enrichment_mc3)

    return x, y
