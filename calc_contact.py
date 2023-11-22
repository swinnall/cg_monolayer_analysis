""" contact map analysis - adapted from Raquel Lopez-Leos De Castro """

import sys
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis import distances
import MDAnalysis as mda
from MDAnalysis.analysis import pca
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA



def contact_matrix(verbose, u, lower_frame, upper_frame, ag1, ag2):

    #ag1: contact group 1 - here you can put your atom selection
    #ag2: contact group 2 - here you can put your atom selection

    # select any resid
    resid_1 = np.unique(ag1.resids)[0]
    resid_2 = np.unique(ag2.resids)[0]

    # select all atoms from input atom group for the above selected resid
    ag1_n = int(len(ag1.select_atoms('resid '+str(resid_1))))
    ag2_n = int(len(ag2.select_atoms('resid '+str(resid_2))))

    # CODE UPDATE - OVERRIDING ag2_n
    # ag2_n = len(ag2)

    contacts_mat= np.full((ag1_n, ag2_n), 0)

    if verbose == True:
        print('\nfirst resid from atom group:',resid_1,resid_2)
        print('number of atoms in one resid:',ag1_n,ag2_n)
        print('contact matrix shape:',np.shape(contacts_mat),'\n')


    for ts in u.trajectory[lower_frame:upper_frame]:
        print("frame=%s/%s" %(ts.frame,len(u.trajectory)-1))

        big_mat_t = mda.analysis.distances.distance_array(ag1.positions,ag2.positions)

        contact_distance = 5 # contact distance

        indices_contact = np.where(big_mat_t <= contact_distance)
        indices_non_contact = np.where(big_mat_t > contact_distance)
        big_mat_t[indices_contact] = 1
        big_mat_t[indices_non_contact] = 0
        big_mat = big_mat_t

        for i in range(0, len(big_mat_t), ag1_n):
            for j in range(0, len(big_mat_t[0]), ag2_n):

                if i != j:

                    c = big_mat[i:i+ag1_n]

                    mat = [row[j:j+ag2_n] for row in c] # selecting minors

                    contacts_mat = contacts_mat + mat

    total_contacts = np.sum(contacts_mat)
    return contacts_mat




def main(u,sysComp,frame0):

    """
    input atom selection for polyA has to be one atom at a time otherwise too much data
    it iterates through via: for each molecule iterate through each base
    so the first 9-10 x values correspond to one polyA molecule (9-10 bases)
    it's 9-10 because P atoms only appear in 9 out of 10 of the bases
    """

    # select whether to print atom group information
    verbose = True
    if verbose == True: print('\n---\nverbose = True.\n')

    # define frame range
    lower_frame = frame0
    upper_frame = 5000

    # string of atoms for atom groups

    # ATOMS DATABASE
    # MC3 HEAD: N1 CN GLA CX" "resname MC3H"
    # POLYA: "BB BB1 BB2 SC1 SC2 SC3 SC4 SC5 SC6" "resname A"

    ag1_resname_str = "resname A"
    ag2_resname_str = "resname A"

    ag1_name_str = "BB BB1 BB2 SC1 SC2 SC3 SC4 SC5 SC6"
    ag2_name_str = "BB BB1 BB2 SC1 SC2 SC3 SC4 SC5 SC6"

    # define strings for atom group selection
    ag1_str = ag1_resname_str + ' and name ' + ag1_name_str
    ag2_str = ag2_resname_str + ' and name ' + ag2_name_str

    # make the atom groups
    ag1 = u.select_atoms(ag1_str)
    ag2 = u.select_atoms(ag2_str)


    if verbose == True:
        print(f'molecule names: {ag1_resname_str}, {ag2_resname_str}')
        print('number of listed atoms:',len(ag1_name_str.split()),len(ag2_name_str.split()))
        print('number of atoms in group:',len(ag1),len(ag2))
        print('number of atom occurences:',len(ag1)/len(ag1_name_str.split()),len(ag2)/len(ag2_name_str.split()))


    # generate contact map
    contacts_mat = contact_matrix(verbose, u, lower_frame, upper_frame, ag1, ag2)


    # define axis labels for heatmap
    y_axis_labels = ag1_name_str.split(' ')
    x_axis_labels = ag2_name_str.split(' ')

    #print(ag2.select_atoms('resid 1').names)


    # index_0 = ag2.indices[0]
    # nAtoms_polyA = len(ag2)
    # x_axis_labels = []
    # for name, resid, index in zip(ag2.names, ag2.resids, ag2.indices):
    #     molecule_id = int( (index - index_0) / nAtoms_polyA)
    #     x_axis_labels.append(f'm{molecule_id}_b{resid}')
    #     # print(f'{molecule_id}_{resid}')
    #
    # x_axis_labels = [i for i in range(len(ag2))]


    if verbose == True:
        print(f'\nnumber of axis labels (x,y): {len(x_axis_labels)}, {len(y_axis_labels)}')


    # fig, ax = plt.subplots(figsize=(15,25)) #you will need to adjust this depending on the shape of your system

    # # number of contacts < 5 A. light color means very few
    # ax = sns.heatmap(contacts_mat, xticklabels=x_axis_labels, yticklabels=y_axis_labels, cmap="YlGnBu")

    # plt.xlabel('PolyA')
    # plt.ylabel('MC3')
    # plt.title('MC3_head-PolyA contacts')

    # plt.show()
    print(np.shape(contacts_mat),contacts_mat,x_axis_labels,y_axis_labels)
    return contacts_mat, x_axis_labels, y_axis_labels



if __name__ == '__main__':
    main()
