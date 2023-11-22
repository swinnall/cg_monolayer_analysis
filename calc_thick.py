import warnings
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import lipyphilic
from lipyphilic.lib.z_thickness import ZThickness


def getThickStatistics(sysComp,membrane,z_thickness):

    cil_str = sysComp[0]

    if 'CHOL' not in sysComp:

        mc3_mask = membrane.resnames == cil_str
        mc3_z_thickness = z_thickness.z_thickness[mc3_mask]

        flat_lipidthick_df_mc3 = pd.DataFrame(mc3_z_thickness[:,:].flatten(), columns=[cil_str])
        print("\n\n%s" %flat_lipidthick_df_mc3.describe())
        print("\nStandard Error:\n%s" %flat_lipidthick_df_mc3.sem())

        mean = [flat_lipidthick_df_mc3.describe().loc[['mean']].iloc[0]]
        se = [flat_lipidthick_df_mc3.sem()]

        return flat_lipidthick_df_mc3, mean, se



    else:

        mc3_mask  = membrane.resnames == cil_str
        chol_mask = membrane.resnames == 'CHOL'

        mc3_z_thickness  = z_thickness.z_thickness[mc3_mask]
        chol_z_thickness = z_thickness.z_thickness[chol_mask]

        flat_lipidthick_df_mc3 = pd.DataFrame(mc3_z_thickness[:,:].flatten(), columns=[cil_str])
        print("\n\n%s" %flat_lipidthick_df_mc3.describe())
        print("\nStandard Error:\n%s" %flat_lipidthick_df_mc3.sem())

        flat_lipidthick_df_chol = pd.DataFrame(chol_z_thickness[:,:].flatten(), columns=['CHOL'])
        print("\n\n%s" %flat_lipidthick_df_chol.describe())
        print("\nStandard Error:\n%s" %flat_lipidthick_df_chol.sem())

        numbers = flat_lipidthick_df_chol['CHOL']
        df = flat_lipidthick_df_mc3.join(numbers)

        mean_mc3  = flat_lipidthick_df_mc3.describe().loc[['mean']].iloc[0]
        mean_chol = flat_lipidthick_df_chol.describe().loc[['mean']].iloc[0]
        mean = [mean_mc3,mean_chol]

        se_mc3  = flat_lipidthick_df_mc3.sem()
        se_chol = flat_lipidthick_df_chol.sem()
        se = [se_mc3,se_chol]

        return df, mean, se



def main(u,sysComp,frame0):

    # filter warnings
    warnings.filterwarnings("ignore")

    # state system
    print('\n\nSystem Components: %s' %(sysComp))

    # define selection criteria
    if 'CHOL' in sysComp:
        lipid_sel_str = 'resname MC3H CHOL'

    else:
        lipid_sel_str = 'resname MC3H'


    # calculate the Z thicknesses
    z_thickness_sn1 = ZThickness(
      universe=u,
      lipid_sel=lipid_sel_str
    )
    z_thickness_sn1.run(start=frame0, stop=None, step=None, verbose=True)

    z_thickness_sn2 = ZThickness(
      universe=u,
      lipid_sel=lipid_sel_str
    )
    z_thickness_sn2.run(start=frame0, stop=None, step=None, verbose=True)

    # combine the thicknesses of each tail into one object
    z_thickness = ZThickness.average(
        z_thickness_sn1,
        z_thickness_sn2
    )

    # define membrane based on assign leaflets text
    membrane = u.select_atoms('name N1 ROH')

    # get dataframe(s) and statistic(s)
    df, mean, se = getThickStatistics(sysComp,membrane,z_thickness)

    return df, mean, se
