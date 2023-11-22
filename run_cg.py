
# standard
import numpy as np
import pandas as pd
import sys
import os
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# gromacs
import MDAnalysis as mda
import lipyphilic
from lipyphilic.lib.assign_leaflets import AssignLeaflets

# for this program
from get_file import getFile, getNRdata
from plot_cg import GenPlot
import calc_apm
import calc_thick
import calc_scc
import calc_enrich
import calc_sld
import calc_rdf
import calc_contact
import get_system


def get_universe(dir):

    XTC = dir + 'p.xtc'

    GRO = dir + 'p.gro'

    TPR = dir + 'p.tpr'

    if os.path.isfile(TPR) == True and os.path.isfile(XTC) == True:
        u = mda.Universe(TPR, XTC, continuous=True)
        print('\nLoaded universe: %s.' %dir)
        print('Trajectory Size: %d frames.' %len(u.trajectory))

    elif os.path.isfile(GRO) == True and os.path.isfile(XTC) == True:
        u = mda.Universe(GRO, XTC, continuous=True)
        print('\nLoaded universe: %s.' %dir)
        print('Trajectory Size: %d frames.' %len(u.trajectory))

    else:
        u = None
        print(f'\nFile in dir ({dir}) does not exist - skipping.\n')

    return u



def get_fig_pars(analysisType):

    if analysisType == 'apm':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = '',
                        yAxisText = r'$\huge{\alpha\ \rm{(\mathring{A}^{2})}}$',
                        saveName  = 'apm',
                        )

    if analysisType == 'thick':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = '',
                        yAxisText = r'$\huge{\tau \ \rm{(\mathring{A})}}$',
                        saveName  = 'thick',
                        )

    if analysisType == 'scc':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = r'$\huge{C_i}$',
                        yAxisText =  r'$\huge{SCC}$',
                        saveName  = "scc",
                        )

    if analysisType == 'enrich':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = "Frame Number",
                        yAxisText =  "Enrichment Index",
                        saveName  = "enrich",
                        )

    if analysisType == 'sld':
        par = dict(x_axisType = 'linear',
                y_axisType = 'linear',
                xAxisText = r'$\huge{\rm{Z \ \left(\mathring{A}\right) }}$',
                yAxisText = r"$\huge{\rm{SLD \ \left(10^{-6} \mathring{A}^{-2}\right) }}$",
                saveName  = "sld",
                )

    if analysisType == 'refl':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'log',
                    xAxisText = r'$\huge{ \rm{Q \ \left(\mathring{A}^{-1}\right) }}$',
                    yAxisText = r'$\huge{\rm{R}}$',
                    saveName  = "R",
                    )

    if analysisType == 'rdf':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'linear',
                    xAxisText = r'$\huge{r \rm{(\mathring{A})}}$',
                    yAxisText = r'$\huge{g(r)}$',
                    saveName  = "rdf",
                    )

    if analysisType == 'contact':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'linear',
                    xAxisText = "PolyA",
                    yAxisText =  "MC3",
                    saveName  = "contact",
                    )

    return par


def get_chol_perc(sysComp):

    # define percentage of cholesterol
    if 'CHOL' not in sysComp:
        cholPerc = '0'
        leaflet_str = "name N1"

    elif 'CHOL' in sysComp and '25' in sysComp:
        cholPerc = '25'
        leaflet_str = 'name N1 ROH'

    elif 'CHOL' in sysComp and '43' in sysComp:
        cholPerc = '43.5'
        leaflet_str = "name N1 ROH"

    # create corresponding string for labelling
    cholName = cholPerc + '% chol.'

    return cholName, leaflet_str



# get root folder
root = get_system.root


# get the system matrix
sysMat = get_system.sysMat


# define number of rows and columns
nRows = np.shape(sysMat)[0]
nCols = np.shape(sysMat)[1]

# choose which analysis to perform: apm / thick / scc / enrich / sld / rdf / contact
analysisType = 'enrich'

# initialise figure class object
# violin plots have one column only, rest have nCols
fig = GenPlot(rows=nRows,cols=1,subplot_titles='')

# get the parameters to define figure axes
fig_par = get_fig_pars(analysisType)

# define second plot if dealing with SLD generation and get NR data
if analysisType == 'sld':
    fig2 = GenPlot(rows=nRows,cols=1,subplot_titles='')
    fig_par2 = get_fig_pars('refl')
    Q, expNR, expNR_err, labels = getNRdata(nCols)

# set the starting frame
frame0 = -50


# iterate through each system in sysmat and call analysis module
for row in range(nRows):
    for col in range(nCols):

        # access name for given figure element
        sysName = sysMat[row][col]

        # create a list of the system components
        sysComp = sysName.split('_')

        # define the file directory
        dir = '../'+root+'/'+sysName+'/'

        # get universe
        u = get_universe(dir)

        # get the amount of cholesterol in the system
        cholName, leaflet_str = get_chol_perc(sysComp)

        # get leaflets
        if analysisType in ['apm'] and u != None:

            print('\nAssigning leaflets:')

            leaflets = AssignLeaflets(
              universe=u,
              lipid_sel=leaflet_str,
            )
            leaflets.run(start=frame0, stop=None, step=None, verbose=True)


        """
        Figure setup
        """

        # reset column to zero for correct plotting
        if analysisType in ['apm','thick','enrich','sld']:
            if analysisType == 'sld':
                Q_ = Q.get(col)
                expNR_ = expNR.get(col)
                expNR_err_ = expNR_err.get(col)
                label_ = labels.get(col)
                col_orig = col
            col = 0

        # define axis labels
        if row == nRows-1:
            xAxisLabel = fig_par.get('xAxisText')
            show_xticklabels = True
        else:
            xAxisLabel = ""
            show_xticklabels = False

        if col == 0:
            yAxisLabel = fig_par.get('yAxisText')
            show_yticklabels = True
        else:
            yAxisLabel = ""
            show_yticklabels = False

        # populate the subplot with the axis parameters
        fig.update_xaxis(xAxisLabel,show_xticklabels,fig_par.get('x_axisType'),row=row+1,col=col+1)
        fig.update_yaxis(yAxisLabel,show_yticklabels,fig_par.get('y_axisType'),row=row+1,col=col+1)


        """
        Analysis modules
        """

        if analysisType == 'apm':
            if u != None:
                df, mean, se = calc_apm.main(u,leaflets,sysComp,frame0)
            else:
                df = pd.DataFrame(columns=[sysComp[0],'CHOL'])
            fig.plotViolin_halves(df,sysComp,cholName,row=row+1,col=col+1,showlegend=False)

        if analysisType == 'thick':
            df, mean, se = calc_thick.main(u,sysComp,frame0)
            fig.plotViolin_halves(df,sysComp,cholName,row=row+1,col=col+1,showlegend=False)

        if analysisType == 'scc':
            calc_scc.main(u,sysName)

        if analysisType == 'enrich':
            x, y = calc_enrich.main(u,sysComp,frame0)
            speciesList = ["feCHL1","feDLMC3"]
            fig.plotLogLog(x,y,sysName,sysComp,speciesList,row=row+1,col=col+1,showlegend=True)
            outputDir = f'output/{fig_par.get("saveName")}_{sysName}'
            fig.fig.write_image(outputDir+'.png',scale=5)

        if analysisType == 'sld':

            SLD_x, SLD_y, contrast, nFrames, q, R = calc_sld.main(u,sysComp,label_,Q_,frame0)

            fig.plotSLD(SLD_x,SLD_y,contrast,nFrames,col_orig,row=row+1,col=col+1,showlegend=True)

            fig2.update_xaxis(fig_par2.get('xAxisText'),show_xticklabels,fig_par2.get('x_axisType'),row=row+1,col=col+1)
            fig2.update_yaxis(fig_par2.get('yAxisText'),show_yticklabels,fig_par2.get('y_axisType'),row=row+1,col=col+1)
            fig2.plotR(contrast,nFrames,q,R,Q_,expNR_,expNR_err_,label_,col_orig,row=row+1,col=col+1,showlegend=True)

            outputDir = f'output/{fig_par.get("saveName")}_{sysName}'
            fig.fig.write_image(outputDir+'.png',scale=5)
            fig2.fig.write_image(outputDir+'_refl.png',scale=5)


        if analysisType == 'rdf':
            x, y = calc_rdf.main(u,sysComp,frame0)
            speciesList = ["MC3"] if 'CHOL' not in sysComp else ["MC3","CHOL"]
            fig.plotLogLog(x,y,sysName,sysComp,speciesList,row=row+1,col=col+1,showlegend=False)



        if analysisType == 'contact':

            contacts_mat, xticklabels, yticklabels = calc_contact.main(u,sysComp,frame0)

            if row == nRows-1:
                xAxisLabel = "PolyA atoms"# (" + atom_str + ")"
                show_xticklabels = True

            else:
                xAxisLabel = ""
                show_xticklabels = False

            if col == 0:
                yAxisLabel = "PolyA atoms" #"MC3 head group atoms"
                show_yticklabels = True
            else:
                yAxisLabel = ""
                show_yticklabels = False

            # populate the subplot with the axis parameters
            fig.update_xaxis(xAxisLabel,show_xticklabels,fig_par.get('x_axisType'),row=row+1,col=col+1)
            fig.update_yaxis(yAxisLabel,show_yticklabels,fig_par.get('y_axisType'),row=row+1,col=col+1)

            fig.plotHeatMap(xticklabels,yticklabels,contacts_mat,row=row+1,col=col+1)

            contact_type = 'polya_polya_mc3_0chol_34polya_10nt'


"""
Show and save figures
"""
if analysisType not in ['scc','enrich','sld','contact']:
    outputDir = f'output/{fig_par.get("saveName")}/{fig_par.get("saveName")}_{sysName}'
    fig.fig.write_image(outputDir+'.png',scale=5)

if analysisType == 'contact':
    outputDir = f'output/{fig_par.get("saveName")}/{fig_par.get("saveName")}_{contact_type}.png'
    fig.fig.write_image(outputDir,scale=5)
