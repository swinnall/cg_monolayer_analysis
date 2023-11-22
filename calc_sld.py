from refnx.reflect import reflectivity
import numpy as np
import sys

"""
To do:

    * overplot experimental SLD with MD SLD
        this would require fitting the experimental data...

    * average both monolayers together

    * generate volume fractions from SLD

"""


class SLD:

    def __init__(self,u,thickness,region_of_interest,z_pos_lower):

        # define scattering lengths (fm)
        self.SL_atoms_AA = {
               "H": -3.739,
               "D": 6.67,
               "C": 6.646,
               "N": 9.36,
               "O": 5.803,
               "P": 5.10,
               }

        # scattering lengths for coarse grain beads [fm]
        self.SL_atoms_CG = {
                # mc3
                    'N1': 0.218,
                    'CN': -1.664,
                    'GLA': 17.42,
                    'CX': 1.243,
                    'C1A': -3.328,
                    'C1Ad': 79.944,
                    'D2A': 0.411,
                    'D2Ad': 73.274,
                    'D3A': 7.889,
                    'D3Ad': 59.934,
                    'C4A': -7.899,
                    'C4Ad': 106.6,
                    'C1B': -3.328,
                    'C1Bd': 79.944,
                    'D2B': 0.411,
                    'D2Bd': 73.274,
                    'D3B': 7.889,
                    'D3Bd': 59.934,
                    'C4B': -7.899,
                    'C4Bd': 106.6,
                # water
                    'W': -6.699,
                    'D': 76.572,
                # chol
                    'ROH':  3.307,
                    'ROHd': 57.757,
                    'R1':   11.628,
                    'R1d':  47.928,
                    'R2':   8.721,
                    'R2d':  35.946,
                    'R3':   4.15,
                    'R3d':  58.6,
                    'R4':   4.982,
                    'R4d':  41.282,
                    'R5':   -4.571,
                    'R5d':  22.654,
                    'R6':   -4.571,
                    'R6d':  22.654,
                    'C1':   -2.496,
                    'C1d':  51.954,
                    'C2':   -7.899,
                    'C2d':  91.926,
                # polya
                    'BB': 28.312,
                    'BB1': 14.524,
                    'BB2': 17.42,
                    'SC1': 8.528,
                    'SC2': 8.528,
                    'SC3': 8.528,
                    'SC4': 1.882,
                    'SC5': 8.528,
                    'SC6': 8.528,
                    }

        # add universe to memory
        self.u = u

        # universe box dimensions
        self.x_dimension, self.y_dimension, self.z_dimension = \
            self.u.dimensions[0], self.u.dimensions[1], self.u.dimensions[2]

        # size of layer slice (slab thickness)
        self.thickness = thickness

        # number of slices
        self.slice_num = int(self.z_dimension/self.thickness)

        # import region of interest
        self.region_of_interest = region_of_interest

        # import lower z position
        self.z_pos_lower = z_pos_lower



    def calcSLD_oneFrame(self,frame,contrast):

        # split contrast information into list of strings
        contrast_components = contrast.split('_')

        # create dict of slices each with list of SLs in that slice
        SL_dict = { i: [] for i in range(self.slice_num)  }

        # set universe to specific frame in the trajectory
        self.u.trajectory[frame]

        # list of atoms to be kept as hydrogenous in the mc3 lipid
        mc3_H_atoms = 'HN11 HN12 HN13 HN21 HN22 HN23 H11 H12 H21 H22 H31 H32 H5\
                        H11R H11S H12R H12S'

        # for charged mc3
        mc3_H_atoms = 'HN4 HN11 HN12 HN13 HN21 HN22 HN23 H11 H12 H21 H22 H31 H32 H5\
                        H11R H11S H12R H12S'

        # for polyA
        polyA_H_atoms = "H61 H62"

        # convert to list
        mc3_H_atoms_list = mc3_H_atoms.split()
        polyA_H_atoms_list = polyA_H_atoms.split()


        # iterate over each atom in universe and add SL to given slice
        count = 0
        for atom in self.u.universe.atoms:

            # get z position of given atom
            z_pos = atom.position[2]

            slice_key = int(z_pos/self.thickness)

            # get element of given atom
            ele = atom.name[0] #if atom.name != "CLA" else "N/A"

            # get scattering length of element
            #SL = self.SL_atoms_AA.get(ele)
            SL = self.SL_atoms_CG.get(atom.name)

            # set SL to 0 if not in dict; e.g. ions
            if SL == None:
                SL = 0

            # get resname of the atom to identify which component is contributing
            resname = atom.resname

            # if atom.name == 'CN': print(resname)

            if resname in ['DLMC3','MC3H', 'MC3'] and slice_key < self.slice_num:

                # only deuterate H atoms
                if 'dlipid' in contrast_components and ele == 'H'\
                    and atom.name not in mc3_H_atoms_list:

                        # swap H for D, assuming 100 % deuteration
                        SL = self.SL_atoms.get("D")

                elif 'dlipid' in contrast_components and \
                    atom.name+'d' in self.SL_atoms_CG:
                    SL = self.SL_atoms_CG.get(atom.name+'d')



                # # only deuterate H atoms
                # if 'dlipid' in contrast_components and ele == 'C'\
                #     and atom.name not in mc3_H_atoms_list:

                #         # swap H for D, assuming 100 % deuteration
                #         SL =




            if resname in ['CHL1','CHOL'] and slice_key < self.slice_num:

                if 'dlipid' in contrast_components and ele == 'H' \
                   and atom.name != 'H3':

                        # swap H for D, multiply by 0.8 for % deuteration factor
                        SL = self.SL_atoms.get("D") * 0.8

                if 'dlipid' in contrast_components and \
                    atom.name+'d' in self.SL_atoms_CG:
                    SL = self.SL_atoms_CG.get(atom.name+'d')
                    # deuteration factor of 0.8 already accounted for in dict



            if resname in ['TIP3','W'] and slice_key < self.slice_num:
                if 'h2o' in contrast_components:
                    pass # this is the default SL configuration

                elif 'acmw' in contrast_components:
                    SL = 0 # acmw has scattering length of 0

                elif 'd2o' in contrast_components and ele == 'H':
                    SL = self.SL_atoms.get("D") # d2o swaps H for D

                elif 'd2o' in contrast_components and ele == 'W':
                    SL = self.SL_atoms_CG.get("D")



            if resname in ['ADE','A'] and slice_key < self.slice_num:

                if 'd2o' in contrast_components and ele == 'H' \
                    and atom.name in polyA_H_atoms_list:

                        # following paper says two exchangeable H atoms:
                        # Zhou; Nucleic Acids Res. 2005; 33(19): 6361-6371
                        # coarse graining approach to determine NA structures from SANS in sol.
                        SL = self.SL_atoms.get("D")

                elif 'd2o' in contrast_components and atom.name == '??' \
                    and atom.name in polyA_H_atoms_list:
                        SL = self.SL_atoms_CG.get(atom.name+'d')


            # always append to the SL_dict dict (as long as in range)
            if slice_key < self.slice_num and slice_key > 0:
                SL_dict[slice_key].append(SL)



        """
        Calculate SLD
        """

        # define SLD dictionary
        SLD_dict = {}

        # divide SL by slice volume, multiply by 10 for units
        for key, value in SL_dict.items():
            SLD_dict[key] = 10 * sum(value)/(self.x_dimension*self.y_dimension*self.thickness)


        """
        Clean data
        """

        # calculate upper z position
        z_pos_upper = self.z_pos_lower + self.region_of_interest

        # get key
        key_lower = int(self.z_pos_lower/self.thickness)
        key_upper = int(z_pos_upper/self.thickness)

        # initialise lists for final storage
        SLD_x, SLD_y = [], []

        # isolate first monolayer only
        idx = 1
        for key, value in list(SLD_dict.items()):
            if key > key_lower and key < key_upper:
                SLD_x.append(idx * self.thickness)
                SLD_y.append(value)
                idx += 1


        return SLD_x, SLD_y



def main(u,sysComp,label_exp,Q_exp,frame0):

    # size of layer slice (slab thickness) Angstrom
    thickness = 5

    # number of frames to average across
    frame1 = frame0 + 10

    # define number of frames for plotting label
    nFrames = frame1 - frame0

    # define z position for region of interest in angstroms
    region_of_interest = 95 # default for acmw 65
    z_pos_lower = 50 # monolayers = 20; polya = 50  / 90


    """
    make H-D substitutions during the analysis based on received experimental file
    strings are constructed from the LABEL of experimental data
    """

    # list of contrasts to be generated
    components = label_exp.upper().split()

    # initialise contrasts string
    contrast = 'dlipid' if 'D-LIPID' in components else 'hlipid'

    # add subphase contrast
    contrast = contrast + '_d2o' if 'D<SUB>2</SUB>O' in components else contrast + '_acmw'

    # add polyA if necessary
    if 'POLYA' in components:
        contrast = contrast + '_polyA'

    # print contrast
    #print('components: %s' %components)
    print('contrast: %s' %contrast)



    """
    iterate through all frames and store the scattering length density
    """

    # each frame has a list that need sto be stored
    SLD_y_allFrames = []

    # initialise class
    S = SLD(u,thickness,region_of_interest,z_pos_lower)

    for frame in range(frame0,frame1):
        print('frame: %d/%d' %(frame+1-frame0,nFrames))

        # calculate SLD
        SLD_x, SLD_y = S.calcSLD_oneFrame(frame,contrast)

        # append totals to list
    SLD_y_allFrames.append(SLD_y)



    """
    average all SLD_y data for a given frame
    """

    # define dict for averaged data, length based on SLD_x
    SLD_y_average = []

    # iterate throughlist of totals
    for i in range(len(SLD_x)):

        val_to_be_averaged = []
        for SLD_list in SLD_y_allFrames:
            val_to_be_averaged.append(SLD_list[i])

        SLD_y_average.append(np.average(val_to_be_averaged))




    """
    generate reflectivity from MD-SLD with refnx
    """

    # changing to refnx notation for clarity
    structure = SLD_y_average

    # define roughness of each layer
    roughness = 0.1

    # define imaginary SLD component
    SLD_im = 0

    # define background based on contrast input
    bkg_ = 9e-6 if 'ACMW' in components else 5e-6

    # last value in experimental Q range
    # q1 = Q_exp.iloc[-1]
    q1 = 0.35

    # define x (q) and y (R) variables
    q = np.linspace(0.01,q1,int(1e3))

    # define empty array for populating with individual slab layers
    slabs = []

    # generate slabs for each structure
    for SLD_re in structure:
        slabs.append([S.thickness, SLD_re, SLD_im, roughness])

    # convert to array for refnx processing
    slabs = np.array(slabs)

    # calculate and store reflectivity
    R = list( reflectivity(q,slabs,bkg=bkg_) )


    return SLD_x, SLD_y_average, contrast, nFrames, q, R





if __name__ == '__main__':
    main()
