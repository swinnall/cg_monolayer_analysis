import MDAnalysis as md
import numpy as np
import contacts_data as contacts_data
import sys
import warnings

# Suppress specific UserWarnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.topology.PDBParser")
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.topology.guessers")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="xdrlib")

## ---------------------------------------------------------------------------##
# USER INPUTS

# define file names of the input pdb structures to get the molecule
# pdb_name = 'mc3h_chol_rna.pdb'
pdb_name = 'l5h_chol_rna.pdb'

# output molecule
# resname = 'MC3H'
# resname = 'LI5H'
# resname = 'CHOL'
resname = 'A'

## ---------------------------------------------------------------------------##

print(f'\nSample info: {contacts_data.name}')

# define directory
dir = '../output/contact/visualisation/'

# Load the universe
u = md.Universe(dir+'input/'+pdb_name)

# Add the betavalues to the system to be sure
u.add_TopologyAttr('tempfactors')

# Select atoms
sel_num = u.select_atoms(f'resname {resname}').resids[0]
sel = u.select_atoms(f'resname {resname} and resid {sel_num}')

# access contact values (each data point is the mean of p1,p2,p3)
contacts = np.array(contacts_data.contacts)

# access associated std with contact values (each data point is STD over p1,p2,p3)
std = np.array(contacts_data.std)

# select a second contacts matrix for subtraction
try:
    contacts2 = np.array(contacts_data.contacts2)
except AttributeError:
    # Handle the case where contacts_2 doesn't exist or isn't defined
    contacts2 = None

# if second contacts matrix exists then subtract it from contacts
if contacts2 is not None and np.array_equal(contacts.shape, contacts2.shape):
    print(f'\nSubtracting contacts from: {contacts_data.name2}')
    contacts = contacts - contacts2

# print(f'\ncontacts = {contacts}')

total_contacts = np.sum(contacts)
print(f'\nTotal Contacts = {total_contacts:.3f} \n-----')

# Sum each row and convert to a list
sums_per_row = np.round(contacts.sum(axis=1).tolist(), 3)

# Sum each column and convert to a list
sums_per_column = np.round(contacts.sum(axis=0).tolist(), 3)

# Print the result
print(f'\nSums per Row = {sums_per_row} \n\nSums per Col = {sums_per_column} \n-----')


# Sum each row and convert to a list
std_per_row = std.sum(axis=1).tolist()

# Sum each column and convert to a list
std_per_column = std.sum(axis=0).tolist()

# Print the result
print(f'\nStd per Row = {std_per_row} \n\nStd per Col = {std_per_column} \n-----')


# Check if the number of selected atoms matches the size of the contact array
sel.n_atoms == contacts.shape[0]

# output values to
if resname == 'A':

    print(f'\nThe sum of sums_per_row = {sum(sums_per_row):.3f}')
    print(f'The sum of std_per_row = {sum(std_per_row):.3f} \n-----')

    # Assign normalised contact values to tempfactors
    try: sel.atoms.tempfactors = sums_per_column / np.max(sums_per_column)
    except ValueError:
        print('\nValueError: Input matrix dimensions and output atoms not equal.\n')
        sys.exit()

else:
    # Assign normalised contact values to tempfactors
    try: sel.atoms.tempfactors = sums_per_row / np.max(sums_per_row)
    except ValueError:
        print('\nValueError: Input matrix dimensions and output atoms not equal.\n')
        sys.exit()

# Write the selected atoms to a PDB file
# sel.atoms.write(dir+'output/'+contacts_data.name+'_'+resname+'.pdb')

# then in VMD:
# source /Volumes/SAM_USB/phd/thesis/analysis/md/cg/output/contact/visualisation/viz.tcl
# 0.5 bond radius for RNA
