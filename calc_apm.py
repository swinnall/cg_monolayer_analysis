import warnings
import pandas as pd
import sys

import lipyphilic
from lipyphilic.lib.area_per_lipid import AreaPerLipid

warnings.filterwarnings("ignore", category=DeprecationWarning)



def genstats(membrane,areas,mask_str):

    # create a mask of resnames in the membrane via input string
    mask = membrane.resnames == mask_str

    # apply mask to the areas previously calculated
    areas = areas.areas[mask]

    # create a dataframe
    df = pd.DataFrame(areas[:,:].flatten(), columns=[mask_str])

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
    print('\n\nSystem lipids: %s' %(lipids))

    # define leaflet string which should be head group atoms 
    leaflet_str = "name GLA "

    # add cholesterol bead to list if in simulation
    if 'CHOL' in lipids:
        leaflet_str = leaflet_str + "ROH"


    print('\nAssigning leaflets:')

    leaflets = AssignLeaflets(
      universe=u,
      lipid_sel=leaflet_str,
    )
    leaflets.run(start=frame0, stop=None, step=None, verbose=True)


    # define selection criteria
    lipid_sel_str = 'name GLA ' # NP N1 ROH

    # add cholesterol bead to list if in simulation
    if 'CHOL' in lipids:
        lipid_sel_str = lipid_sel_str + 'ROH'

    # Calculate the area per lipid in each frame
    areas = AreaPerLipid(
      universe=u,
      lipid_sel=lipid_sel_str,
      leaflets=leaflets.leaflets
    )
    areas.run(start=frame0, stop=None, step=None, verbose=True)

    # define membrane based on assign leaflets text
    membrane = u.select_atoms(lipid_sel_str)

    # initialise dataframe
    df = pd.DataFrame()

    # for each lipid component in the system, generate statistics and store values
    for lipid in lipids:
            df[lipid] = genstats(membrane,areas,lipid)

    return df
