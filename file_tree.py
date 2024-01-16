"""
list of files for a given input directory

choose systems to studyy:
  [ [row0col0, row0col1 ...], [row1col0, row1col1 ...]  ]
"""

# SLD:
# calcSLD.py will generate the correct contrast based on sample-info-NR.txt
# from the a single MD universe which must be called as many times as there
# are experimental NR contrasts


## ---------------------------------------------------------------------------##
# dir = 'mc3_mono/0p'

# FILES:
# mc3_chol_57_43
# mc3_chol_75_25
# mc3_100

# dir = 'mc3_mono/10p'

# FILES:
# mc3h_mc3_chol_6_51_43
# mc3h_mc3_chol_7.5_67.5_25
# mc3h_mc3_10_90

# dir = 'mc3_mono/100p'

# FILES:
# mc3h_chol_57_43 # one frame??
# mc3h_chol_75_25
# mc3h_100

# mat = [['mc3h_chol_57_43','mc3h_chol_57_43','mc3h_chol_57_43']]


## ---------------------------------------------------------------------------##
# dir = 'mc3_rna_10nt/10p'

# FILES:
# mc3h_mc3_chol_6_51_43_rna_34
# mc3h_mc3_chol_7.5_67.5_25_rna_34
# mc3h_mc3_chol_10_90_0_rna_34


# dir = 'mc3_rna_10nt/34mol'

# FILES:
# mc3_chol_57_43_rna_34
# mc3_chol_75_25_rna_34
# mc3_chol_100_0_rna_34
# mc3h_chol_57_43_rna_34
# mc3h_chol_75_25_rna_34
# mc3h_chol_100_0_rna_34


# dir = 'mc3_rna_10nt/338mol'

# FILES:
# MC3H_CHOL_57_43_RNA_338
# MC3H_CHOL_75_25_RNA_338
# MC3H_CHOL_100_0_RNA_338

# mat = [['mc3h_chol_57_43_rna_34','mc3h_chol_57_43_rna_34','mc3h_chol_57_43_rna_34']]


## ---------------------------------------------------------------------------##
# dir = 'mc3_rna_20nt/100p'

# FILES:
# mc3_chol_57_43_rna_34
# mc3_chol_75_25_rna_34
# mc3_chol_100_0_rna_34
# MC3H_CHOL_57_43_34RNA
# MC3H_CHOL_75_25_34RNA
# MC3H_CHOL_100_0_34RNA

# mat = [['MC3H_CHOL_57_43_34RNA','MC3H_CHOL_57_43_34RNA','MC3H_CHOL_57_43_34RNA']]

## ---------------------------------------------------------------------------##
# dir = 'lipid5_mono/100p'

# FILES:
# LI5H_CHOL_57_43
# LI5H_CHOL_75_25
# LI5H_CHOL_100_0
#
# mat = [['LI5H_CHOL_57_43','LI5H_CHOL_57_43','LI5H_CHOL_57_43']]

dir = 'lipid5_rna_10nt/100p'

# FILES:
# LI5H_CHOL_57_43_RNA
# LI5H_CHOL_75_25_RNA
# LI5H_CHOL_100_0_RNA

mat = [['LI5H_CHOL_57_43_RNA','LI5H_CHOL_57_43_RNA','LI5H_CHOL_57_43_RNA']]
