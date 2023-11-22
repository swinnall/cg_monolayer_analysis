## SYSTEMS TO STUDY

# MC3_monolayers

# MC3
# MC3_CHOL_25
# MC3_CHOL_43

# MC3H
# MC3H_CHOL_25
# MC3H_CHOL_43


# MC3_monolayers_RNA_10nt (34 RNA molecules unless otherwise stated)

# MC3_RNA
# MC3_CHOL_25_RNA
# MC3_CHOL_43_RNA

# MC3H_10P_RNA
# MC3H_10P_CHOL_25_RNA
# MC3H_10P_CHOL_43_RNA

# MC3H_RNA
# MC3H_CHOL_25_RNA
# MC3H_CHOL_43_RNA

# MC3H_RNA_338


# MC3_monolayers_RNA_20nt (34 RNA molecules)

# MC3H_34RNA
# MC3H_CHOL_25_34RNA
# MC3H_CHOL_43_34RNA


# ---------------------------------------------------------------------------- #
# choose root folder
root = 'MC3_monolayers'

# choose systems to studyy: [ [row0col0, row0col1 ...], [row1col0, row1col1 ...]  ]

# sysMat = [['MC3_RNA']]

# sysMat = [['MC3H_RNA_338']]

# sysMat = [['MC3H','MC3H_CHOL_25','MC3H_CHOL_43']]
sysMat = [['MC3H_CHOL_25','MC3H_CHOL_43']]

# sysMat = [['MC3H_RNA','MC3H_CHOL_25_RNA','MC3H_CHOL_43_RNA']]
# sysMat = [['MC3H_CHOL_25_RNA','MC3H_CHOL_43_RNA']]

# sysMat = [['MC3H_34RNA','MC3H_CHOL_25_34RNA','MC3H_CHOL_43_34RNA']]
# sysMat = [['MC3H_CHOL_25_34RNA','MC3H_CHOL_43_34RNA']]

# sysMat = [['MC3H_CHOL_25_RNA']]

# sysMat = [['MC3','MC3_CHOL_25','MC3_CHOL_43']]

# sysMat = [['MC3','MC3_CHOL_25','MC3_CHOL_43'],\
          # ['MC3H','MC3H_CHOL_25','MC3H_CHOL_43']]

# sysMat = [['MC3H_10P_RNA','MC3H_10P_CHOL_25_RNA']]
# sysMat = [['MC3H_34RNA','MC3H_CHOL_25_34RNA','MC3H_CHOL_43_34RNA']]

# sysMat = [['MC3_RNA','MC3_CHOL_25_RNA','MC3_CHOL_43_RNA'],\
#           ['MC3H_10P_RNA','MC3H_10P_CHOL_25_RNA','MC3H_10P_CHOL_43_RNA'],\
#           ['MC3H_RNA','MC3H_CHOL_25_RNA','MC3H_CHOL_43_RNA']]






# ---------------------------------------------------------------------------- #
# SLD:
# calcSLD.py will generate the correct contrast based on sample-info-NR.txt
# from the a single MD universe which must be called as many times as there
# are experimental NR contrasts

# sysMat = [['MC3_CHOL_25_RNA','MC3_CHOL_25_RNA','MC3_CHOL_25_RNA']]

# sysMat = [['MC3H_RNA_338','MC3H_RNA_338','MC3H_RNA_338']]

#sysMat = [['cg_mc3_chol25','cg_mc3_chol25','cg_mc3_chol25']]
#sysMat = [['cg_mc3_chol43','cg_mc3_chol43','cg_mc3_chol43']]
