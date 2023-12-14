
from lipyphilic.lib.order_parameter import SCC
from lipyphilic.lib.assign_leaflets import AssignLeaflets
from lipyphilic.lib.plotting import ProjectionPlot
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# SCD is the order parameter for atomistic simulations
# SCC is for coarse grain and is the secondary order parameter

def main(u,fname,frame0,rep):

    sel_str = 'resname MC3H MC3' # 'resname MC3H CHOL' name N1 NP

    print('\nAssigning leaflets:')
    leaflets = AssignLeaflets(universe=u,lipid_sel=sel_str)
    leaflets.run(start=frame0,stop=None,step=None,verbose=True)

    print('\nCalculating tail 1')
    scc_sn1 = SCC(universe=u,tail_sel="name ??A") # selects C1A, C2A, D2A, C3A, and C4A
    scc_sn1.run(start=frame0,stop=None,step=1,verbose=True)

    print('\nCalculating tail 2')
    scc_sn2 = SCC(universe=u,tail_sel="name ??B")  # selects C1B, C2B, D2B, C3B, and C4B
    scc_sn2.run(start=frame0,stop=None,step=1,verbose=True)

    # returns np array, row = individual lipid, col = individual frame
    scc_av = SCC.weighted_average(scc_sn1, scc_sn2)

    # print range of data
    print(f'\nMin = {np.min(scc_av.SCC)}\nMax = {np.max(scc_av.SCC)}')

    # average across the lipids (averaging column rows), enables plotting scc against time
    scc_av_lipid = np.mean(scc_av.SCC,axis=0)

    # average across the frames average array
    scc_av_lipid_frame = np.mean(scc_av_lipid)

    # calculate the standard deviation
    scc_std_lipid_frame = np.std(scc_av_lipid)

    print(f'\nAverage SCC = {scc_av_lipid_frame:.3f}')
    print(f'Std SCC = {scc_std_lipid_frame:.3f}\n-----\n')

    for leaflet_idx in range(2):

        # cbar = False if leaflet_idx == 0 else True
        cbar = True

        leaflet_num = 1 if leaflet_idx == 0 else -1

        leaflet_str = 'upper' if leaflet_idx == 0 else 'lower'

        scc_projection = scc_av.project_SCC(
              lipid_sel=sel_str,
              start=frame0,
              stop=None,
              step=None,
              filter_by=leaflets.filter_leaflets('name ??A ??B') == leaflet_num,
              cmap='viridis',
              cbar=cbar,
              vmin=-0.05,
              vmax=0.1,
              bins=100,
              cbar_kws=dict(shrink=0.83)
        )

        scc_projection.fig.savefig(f'../output/scc/scc_p{rep}_{fname}_{leaflet_str}.png',dpi=300)

    return
