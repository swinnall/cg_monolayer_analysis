
# standard
import numpy as np
import pandas as pd
import sys
import os
import warnings
import MDAnalysis as mda

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
import file_tree


def get_universe(dir,fname,rep):

    path = dir + '/' + fname + '/'

    stride = False
    for fname in os.listdir(path):
        if 'stride' in fname:
            stride = True


    # define file names
    if rep != None:

        XTC = path + f'p_{rep}_stride.xtc' if stride == True else path + f'p_{rep}.xtc'
        GRO = path + f'p_{rep}.gro'
        TPR = path + f'p_{rep}.tpr'

    else:
        XTC = path + 'p.xtc'
        GRO = path + 'p.gro'
        TPR = path + 'p.tpr'

    # load universe
    if os.path.isfile(TPR) == True and os.path.isfile(XTC) == True:
        u = mda.Universe(TPR, XTC, continuous=True)
        print('\nLoaded MD files: \n %s \n %s' %(TPR,XTC))
        print('Trajectory Size: %d frames.' %len(u.trajectory))

    elif os.path.isfile(GRO) == True and os.path.isfile(XTC) == True:
        u = mda.Universe(GRO, XTC, continuous=True)
        print('\nLoaded MD files: \n %s \n %s' %(GRO,XTC))
        print('Trajectory Size: %d frames.' %len(u.trajectory))

    else:
        u = None
        print(f'\nFile in dir ({path}) does not exist - skipping.\n')

    return u



def get_fig_pars(analysis):

    if analysis == 'apm':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = '',
                        yAxisText = r'$\huge{\alpha\ \rm{(\mathring{A}^{2})}}$',
                        saveName  = 'apm',
                        )

    if analysis == 'thick':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = '',
                        yAxisText = r'$\huge{\tau \ \rm{(\mathring{A})}}$',
                        saveName  = 'thick',
                        )

    if analysis == 'scc':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = r'$\huge{C_i}$',
                        yAxisText =  r'$\huge{SCC}$',
                        saveName  = "scc",
                        )

    if analysis == 'enrich':
        par = dict(x_axisType = 'linear',
                        y_axisType = 'linear',
                        xAxisText = r"$\huge{\text{Frame Number}}$",
                        yAxisText =  r"$\huge{\text{Enrichment Index}}$",
                        saveName  = "enrich",
                        )

    if analysis == 'sld':
        par = dict(x_axisType = 'linear',
                y_axisType = 'linear',
                xAxisText = r'$\huge{\rm{Z \ \left(\mathring{A}\right) }}$',
                yAxisText = r"$\huge{\rm{SLD \ \left(10^{-6} \mathring{A}^{-2}\right) }}$",
                saveName  = "sld",
                )

    if analysis == 'refl':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'log',
                    xAxisText = r'$\huge{ \rm{Q \ \left(\mathring{A}^{-1}\right) }}$',
                    yAxisText = r'$\huge{\rm{R}}$',
                    saveName  = "R",
                    )

    if analysis == 'rdf':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'linear',
                    xAxisText = r'$\huge{r \rm{(\mathring{A})}}$',
                    yAxisText = r'$\huge{g(r)}$',
                    saveName  = "rdf",
                    )

    if analysis == 'contact':
        par = dict(x_axisType = 'linear',
                    y_axisType = 'linear',
                    xAxisText = "PolyA",
                    yAxisText =  "MC3",
                    saveName  = "contact",
                    )

    return par


## ---------------------------------------------------------------------------##

# choose which analysis to perform: apm / thick / scc / enrich / sld / rdf / contact
analysis = 'contact'

# set the starting frame
frame0 = -25 #-25

## ---------------------------------------------------------------------------##


# disable warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# keep these words if in components
keep_lipid = set(['MC3','MC3H','LIPID5','LIPID5H','LI5H','CHOL'])

# keep these words if in components
keep_amount = set(['25','43'])

# get root folder
dir = '../input/' + file_tree.dir

# get the system matrix
fmat = file_tree.mat

# define number of rows and columns
nRows = np.shape(fmat)[0]
nCols = np.shape(fmat)[1]

# initialise figure class object
# violin plots have one column only, rest have nCols
fig = GenPlot(rows=nRows,cols=1,subplot_titles='')

# get the parameters to define figure axes
fig_par = get_fig_pars(analysis)

# define second plot if dealing with SLD generation and get NR data
if analysis == 'sld':
    fig2 = GenPlot(rows=nRows,cols=1,subplot_titles='')
    fig_par2 = get_fig_pars('refl')
    Q, expNR, expNR_err, labels = getNRdata(nCols)

list_enrich, list_contacts_mat = [], []

# iterate through each system in sysmat and call analysis module
for row in range(nRows):
    for col in range(nCols):

        # access name for given figure element
        fname = fmat[row][col]

        # create a list of the system components
        components = fname.split('_')

        # make everything capital letters
        components = [str.upper() for str in components]

        # make list of lipids from the input components
        lipids = [str for str in components if str in keep_lipid]

        # extract percentage cholesterol in simulation, should be one value
        chol_percentage = [num for num in components if num in keep_amount]

        # define cholesterol name
        if chol_percentage != []:
            cholName = chol_percentage[0] + '% chol.'
        else:
            cholName = '0% chol.'

        # define rep if there are repeats or not: row+1 / None
        rep = col+1

        # ensure using the same (p1) input file for the sld calculations
        if analysis == 'sld': rep = 2

        # get universe
        u = get_universe(dir,fname,rep)


        """
        Figure setup
        """

        # reset column to zero for correct plotting
        if analysis in ['apm','thick','enrich','sld']:
            if analysis == 'sld':
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

        if analysis == 'apm':
            if u != None:
                df = calc_apm.main(u,lipids,frame0)
            else:
                df = pd.DataFrame(columns=[lipids[0],'CHOL'])
            fig.plotViolin(df,lipids,cholName,row=row+1,col=col+1,showlegend=False)

        if analysis == 'thick':
            if u != None:
                df = calc_thick.main(u,lipids,frame0)
            else:
                df = pd.DataFrame(columns=[lipids[0],'CHOL'])
            fig.plotViolin(df,lipids,cholName,row=row+1,col=col+1,showlegend=False)

        if analysis == 'scc':
            calc_scc.main(u,fname,frame0,rep)

        if analysis == 'enrich':
            df = calc_enrich.main(u,lipids,frame0)
            list_enrich.append(df)

        if analysis == 'sld':

            SLD_x, SLD_y, contrast, nFrames, q, R = calc_sld.main(u,lipids,label_,Q_,frame0)

            fig.plotSLD(SLD_x,SLD_y,contrast,nFrames,col_orig,row=row+1,col=col+1,showlegend=True)

            fig2.update_xaxis(fig_par2.get('xAxisText'),show_xticklabels,fig_par2.get('x_axisType'),row=row+1,col=col+1)
            fig2.update_yaxis(fig_par2.get('yAxisText'),show_yticklabels,fig_par2.get('y_axisType'),row=row+1,col=col+1)
            fig2.plotR(contrast,nFrames,q,R,Q_,expNR_,expNR_err_,label_,col_orig,row=row+1,col=col+1,showlegend=True)

            outputDir = f'../output/{fig_par.get("saveName")}_{fname}'
            fig.fig.write_image(outputDir+'.png',scale=5)
            fig2.fig.write_image(outputDir+'_refl.png',scale=5)


        if analysis == 'rdf':
            x, y = calc_rdf.main(u,lipids,frame0)
            speciesList = ["MC3"] if 'CHOL' not in lipids else ["MC3","CHOL"]
            fig.plotLogLog(x,y,fname,lipids,speciesList,row=row+1,col=col+1,showlegend=False)



        if analysis == 'contact':

            # Select atoms based on the specified resnames
            ag_cil = u.select_atoms(f"resname MC3H LI5H")
            ag_chol = u.select_atoms(f"resname CHOL")
            ag_rna = u.select_atoms(f"resname A")

            # Get the number of different residue IDs
            nCIL = len(set(ag_cil.resids))
            nCHOL = len(set(ag_chol.resids))
            nNT = len(set(ag_rna.resids)) # note each nucleotide has different resid)

            # setting the normalisation factors
            # normFactor = nNT; normString = 'nNT' # RNA - RNA
            # normFactor = nCHOL * nNT; normString = 'nChol * nNT' # CHOL - RNA
            normFactor = nCIL * nNT; normString = 'nCIL * nNT' # CIL - RNA
            # normFactor = nCIL * nCHOL; normString = 'nCIL * nCHOL' # CIL - CHOL

            print(f'\nContact Normalisation\nnCIL={nCIL},nChol={nCHOL},nNucleo={nNT}\nnormalisation={normString} = {normFactor}')

            contacts_mat, xticklabels, yticklabels = calc_contact.main(u,lipids,frame0)

            # store matrix list in a list for later use
            list_contacts_mat.append(contacts_mat)

            if row == nRows-1:
                xAxisLabel = "RNA Beads"
                # xAxisLabel = "Chol Beads"
                show_xticklabels = True

            else:
                xAxisLabel = ""
                show_xticklabels = False

            if col == 0:
                # yAxisLabel = "RNA Beads"
                yAxisLabel = "MC3 Beads"
                # yAxisLabel = "Lipid5 Beads"
                # yAxisLabel = "Chol Beads"
                show_yticklabels = True
            else:
                yAxisLabel = ""
                show_yticklabels = False

            # populate the subplot with the axis parameters
            fig.update_xaxis(xAxisLabel,show_xticklabels,fig_par.get('x_axisType'),row=row+1,col=col+1)
            fig.update_yaxis(yAxisLabel,show_yticklabels,fig_par.get('y_axisType'),row=row+1,col=col+1)



"""
Show and save figures
"""
if analysis not in ['scc','enrich','sld','contact']:
    outputDir = f'../output/{fig_par.get("saveName")}/{fig_par.get("saveName")}_{fname}'
    fig.fig.write_image(outputDir+'.png',scale=5)



if analysis == 'contact':

    # Calculate the average and standard deviation
    Z_av = np.array(np.mean(np.array(list_contacts_mat), axis=0) / normFactor )
    Z_std = np.array(np.std(np.array(list_contacts_mat), axis=0) / normFactor )

    # fig.plotHeatMap(xticklabels,yticklabels,Z_av,Z_std,row=row+1,col=col+1)

    # contact_type = 'chol_rna_mc3h_43chol_34polya_10nt'
    # outputDir = f'../output/{fig_par.get("saveName")}/{fig_par.get("saveName")}_{contact_type}.png'
    # fig.fig.write_image(outputDir,scale=5)

    print('\nav:')
    print(*np.round(Z_av,5))
    print('\nstd:')
    print(*np.round(Z_std,5))




if analysis == 'enrich':

    fig.plotEnrich(list_enrich,fname,lipids,row=row+1,col=col+1,showlegend=True)

    outputDir = f'../output/enrich/{fig_par.get("saveName")}_{fname}'
    fig.fig.write_image(outputDir+'.png',scale=5)
