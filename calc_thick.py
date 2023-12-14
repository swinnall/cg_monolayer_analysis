import warnings
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import lipyphilic
from lipyphilic.lib.z_thickness import ZThickness



def genstats(membrane,z_thickness,mask_str):

    # create a mask of resnames in the membrane via input string
    mask = membrane.resnames == mask_str

    # apply mask to the z_thickness previously calculated
    thick = z_thickness.z_thickness[mask]

    print('\nLipid:',mask_str)

    # create a dataframe
    df = pd.DataFrame(thick[:,:].flatten(), columns=[mask_str])

    # calculate the mean of the dataframe
    mean = df.describe().loc[['mean']].iloc[0]

    # calculate the standard error of the dataframe
    se = df.sem()

    # print full 'box plot' statistics
    print("\n\n%s" %df.describe())

    # print the standard error of the system
    print("\nStandard Error:\n%s" %se)

    return df



def main(u,lipids,frame0):

    # filter warnings
    warnings.filterwarnings("ignore")

    # state system
    print('\n\nSystem Components: %s' %lipids)

    # define selection criteria
    sel_str = "resname MC3 MC3H" # (name ??1 ??A)

    if 'CHOL' in lipids:
        sel_str = sel_str + ' CHOL' # or (resname CHOL and not name ROH)

    # calculate the z thicknesses of tail 1
    print('\nTail 1:')
    z_thickness_sn1 = ZThickness(universe=u,lipid_sel='name NP GLA ??A') #sel_str
    z_thickness_sn1.run(start=frame0,stop=None,step=None,verbose=True)

    # calculate the z thicknesses of tail 2
    print('\nTail 2:')
    z_thickness_sn2 = ZThickness(universe=u,lipid_sel='name NP GLA ??B')
    z_thickness_sn2.run(start=frame0,stop=None,step=None,verbose=True)

    # combine the thicknesses of each tail into one object
    print('\nCalculating average...')
    z_thickness = ZThickness.average(z_thickness_sn1,z_thickness_sn2)
    print('Complete.')

    # define membrane selection
    memb_str = 'name NP' # N1 NP ROH
    if 'CHOL' in lipids:
        memb_str = memb_str + ' ROH' # or (resname CHOL and not name ROH)

    # create membrane atom group
    membrane = u.select_atoms(memb_str)

    # initialise dataframe
    df = pd.DataFrame()

    # for each lipid component in the system, generate statistics and store values
    for lipid in lipids:
            df[lipid] = genstats(membrane,z_thickness,lipid)

    return df
