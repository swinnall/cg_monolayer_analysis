import warnings
import pandas as pd
import sys

import lipyphilic
from lipyphilic.lib.area_per_lipid import AreaPerLipid


def getAreaStatistics(sysComp,membrane,areas):

    cil_str = sysComp[0]

    if 'CHOL' not in sysComp:

        mc3_mask = membrane.resnames == cil_str
        mc3_areas = areas.areas[mc3_mask]

        flat_lipidAPM_df_mc3 = pd.DataFrame(mc3_areas[:,:].flatten(), columns=[cil_str])
        print("\n\n%s" %flat_lipidAPM_df_mc3.describe())
        print("\nStandard Error:\n%s" %flat_lipidAPM_df_mc3.sem())

        mean = [flat_lipidAPM_df_mc3.describe().loc[['mean']].iloc[0]]
        se = [flat_lipidAPM_df_mc3.sem()]

        return flat_lipidAPM_df_mc3, mean, se


    else:

        mc3_mask  = membrane.resnames == cil_str
        chol_mask = membrane.resnames == "CHOL"

        mc3_areas  = areas.areas[mc3_mask]
        chol_areas = areas.areas[chol_mask]

        flat_lipidAPM_df_mc3 = pd.DataFrame(mc3_areas[:,:].flatten(), columns=[cil_str])
        print("\n\n%s" %flat_lipidAPM_df_mc3.describe())
        print("\nStandard Error:\n%s" %flat_lipidAPM_df_mc3.sem())

        flat_lipidAPM_df_chol = pd.DataFrame(chol_areas[:,:].flatten(), columns=['CHOL'])
        print("\n\n%s" %flat_lipidAPM_df_chol.describe())
        print("\nStandard Error:\n%s" %flat_lipidAPM_df_chol.sem())

        numbers = flat_lipidAPM_df_chol['CHOL']
        df = flat_lipidAPM_df_mc3.join(numbers)

        mean_mc3  = flat_lipidAPM_df_mc3.describe().loc[['mean']].iloc[0]
        mean_chol = flat_lipidAPM_df_chol.describe().loc[['mean']].iloc[0]
        mean = [mean_mc3,mean_chol]

        se_mc3  = flat_lipidAPM_df_mc3.sem()
        se_chol = flat_lipidAPM_df_chol.sem()
        se = [se_mc3,se_chol]

        return df, mean, se




def main(u,leaflets,sysComp,frame0):

    # filter warnings
    warnings.filterwarnings("ignore")

    # state system
    print('\n\nSystem Components: %s' %(sysComp))

    # define selection criteria
    if 'CHOL' not in sysComp:
        lipid_sel_str = 'name N1' # and not resname ADE

    else:
        lipid_sel_str = 'name N1 ROH'


    # Calculate the area per lipid in each frame
    areas = AreaPerLipid(
      universe=u,
      lipid_sel=lipid_sel_str,
      leaflets=leaflets.leaflets
    )
    areas.run(start=frame0, stop=None, step=None, verbose=True)


    # define membrane based on assign leaflets text
    membrane = u.select_atoms(lipid_sel_str)


    # get dataframe(s) and statistic(s)
    df, mean, se = getAreaStatistics(sysComp,membrane,areas)


    return df, mean, se


if __name__ == '__main__':
    main()
