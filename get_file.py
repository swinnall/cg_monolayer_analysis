import pandas as pd


def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def getNRdata(nCols):

    """
    Read sample information from instructions file
    """

    instructionsDir = 'nr_sample_info.txt'

    instructionsFile = getFile(path=instructionsDir, nSkip=1, delim=',')

    nFiles = len(instructionsFile)

    file_paths, inputLabels, contrastList = {}, {}, []

    dir = '../../../nr/input/RB2220338-MC3/ASCII_Feb23_sw2/'

    for i in range(nFiles):

        contrast = instructionsFile["contrast"][i]
        contrastList.append(contrast)

        file_paths[contrast] = dir + instructionsFile["fname"][i] + ".txt"
        inputLabels[contrast] = instructionsFile["label"][i]

    print(f'\nLoaded NR files: \n{file_paths}')

    """
    Read experimental data
    """

    # define number of models to fit
    nExpFiles = nCols

    # initialise variables
    Q, expNR, expNR_err, labels = {}, {}, {}, {}

    # for each contrast, store the experimental information
    for i in range(nExpFiles):
        contrast = contrastList[i]

        if file_paths.get(contrast) is not None:
            experimentFile = getFile(path=file_paths.get(contrast), nSkip=1, delim='\t')

            Q[i] = experimentFile[experimentFile.columns.values[0]]
            expNR[i] = experimentFile[experimentFile.columns.values[1]]
            expNR_err[i] = experimentFile[experimentFile.columns.values[2]]
            labels[i] = inputLabels.get(contrast)

    return Q, expNR, expNR_err, labels
