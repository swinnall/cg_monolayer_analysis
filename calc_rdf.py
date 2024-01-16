import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

import MDAnalysis as mda
from MDAnalysis.analysis import rdf


def get_rdf_statistics(u,ag1,ag2,frame0):

    # create the InterRDF object, by supplying two AtomGroups
    # rdf = numpy.ndarray of the radial distribution function values for the bins
    rdf_ = rdf.InterRDF(ag1,
                        ag2,
                        nbins=75,  # default
                        range=(0.0, 15.0),  # 15 is default distance
                        )

    # perform the calculation
    rdf_.run(start=frame0, stop=None, step=None, verbose=True)

    res = rdf_.results
    rdf_res = res.rdf
    bins = res.bins

    ## find first minimum of the rdf (first hydration shell)

    # get idx of rdf.bins
    first_shell_cutoff = 6.9 # A (originally was 3.7 A but that didn't cover CG so made bigger...)

    cutoff_idx = int(np.where(bins == first_shell_cutoff)[0]) + 1

    # find maximum within set region for the first shell
    max_rdf = max(rdf_res[:cutoff_idx])
    print(cutoff_idx,max_rdf)
    max_idx = int(np.where(rdf_res == max_rdf)[0])

    print('max rdf: %f' %max_rdf)
    print('index of max rdf: %s' %max_idx)


    ### override max index ###
    force_max = False
    if force_max == True:
        max_idx = 14
        print('Overriding automatic max. max_idx = %d' %max_idx)


    # find the minimum nearest to that maximum and before the cutoff
    min_rdf = min(rdf_res[max_idx:cutoff_idx])
    min_idx = int(np.where(rdf_res == min_rdf)[0])

    print('min rdf: %f' %min_rdf)
    print('index of min rdf: %s' %min_idx)

    # get the corresponding radius value at this index
    radius_shell_1 = bins[min_idx]
    print('Radius of first hydration shell = %.1f (A)' %radius_shell_1)


    ### override min index ###
    force_min = False
    if force_min == True:
        min_idx = 17
        print('Overriding automatic min. min_idx = %d' %min_idx)


    # get count in each bin
    count_ = res.count

    # sum over all the bins up to the min index + 1 in order to include the final bin
    count_sum = 0
    for bin_idx in range(min_idx+1):
        print('bin radius = %.1f' %bins[bin_idx])
        count_sum += count_[bin_idx]

    #print('\nTotal Count = %d' %count_sum)

    # total number of atoms across both leaflets -- check for each file??
    nAtoms = 1352

    # calculate the AVERAGE number of molecules in the FIRST hydration shell
    nMolecules = count_sum / (len(u.trajectory) * nAtoms)

    return rdf_, nMolecules



def main(u,sysComp,frame0):

    """
    This calculates the average radial distribution function
    As in: https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/rdf.html#module-MDAnalysis.analysis.rdf
    It is calculated by histogramming distances between all particles in ag1 and ag2
    while taking periodic boundary conditions into account via the minimum image convention.
    """

    # initialise variables for x and y for final plot
    x, y = [], []

    # define atom groups
    # N1 and O3 for MC3 and chol headgroups respectively
    # O1 and O2 for the MC3 carboxylic groups
    # OH2 for TIP3 water molecules
    ag_mc3 = u.select_atoms('resname MC3H and name N1')
    ag_chol = u.select_atoms('resname CHOL and name ROH')
    ag_water = u.select_atoms('resname W')

    print('\nComponent: MC3')
    rdf_mc3, nMolecules_mc3 = get_rdf_statistics(u,ag1=ag_mc3,ag2=ag_water,frame0=frame0)
    x.append(rdf_mc3.results.bins)
    y.append(rdf_mc3.results.rdf)
    print('<N_water/mc3> = %f' %nMolecules_mc3)

    # analyse second component if needed
    if 'CHOL' in sysComp:
        print('\nComponent: Cholesterol')

        rdf_chol, nMolecules_chol = get_rdf_statistics(u,ag1=ag_chol,ag2=ag_water,frame0=frame0)
        x.append(rdf_chol.results.bins)
        y.append(rdf_chol.results.rdf)

        print('<N_water/chol> = %f' %nMolecules_chol)


    quickPlot = False
    if quickPlot == True:

        # plot the radial distribution functions
        plt.plot(rdf_mc3.results.bins, rdf_mc3.results.rdf,label='mc3')

        if 'CHOL' in sysComp:
            plt.plot(rdf_chol.results.bins, rdf_chol.results.rdf,label='chol')

        plt.xlabel('Radius (angstrom)')
        plt.ylabel('Radial distribution')
        plt.legend()
        plt.savefig('output/rdf_quickplot.png')

    return x, y



if __name__ == '__main__':
    main()
