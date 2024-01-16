import MDAnalysis as mda
import os
import sys
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# define output file extension
ext = 'pdb'

# define starting directory
dir = '../input/lipid5_rna_10nt/'

# walk through input file directory
for root, dirs, files in os.walk(dir):
    for name in files:

        # filter out only non hidden tpr files
        if name[-3:] == 'tpr' and name[0] != '.':

            # define input file path name
            TPR = os.path.join(root, name)

            # print file input to terminal for reference
            print(f'\nInput = {TPR}')

            # try different names for xtc files
            XTC1 = root + '/' + name[:-3] + 'xtc'
            XTC2 = root + '/' + name[:-4] + '_stride.xtc'

            # check XTC file exists and select if True
            if os.path.isfile(XTC1) == True:
                XTC = XTC1
            elif os.path.isfile(XTC2) == True:
                XTC = XTC2
            else:
                print(f'\nError: No XTC file exists in current dir. \n{os.listdir(dir)}')

            # load in the universe
            try:
                u = mda.Universe(TPR,XTC)

                # create atom group including all atoms
                ag = u.select_atoms("all")

                # define output file name
                name_out = root + '/' + name[:-3] + ext

                # write new file
                ag.write(name_out)

                # print update to terminal
                print(f'File created = {name_out}')

            except ValueError:
                print(f'Exception: ValueError\nNumber of atoms may not be equivalent in TPR and XTC.\nSkipping...')
