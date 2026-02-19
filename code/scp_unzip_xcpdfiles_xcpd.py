
import os
import re
import glob
import zipfile
import numpy as np
# import datalad.api as dl

dataset = 'HCPYA'

if dataset == 'HCPD':
    datapath = '/cbica/projects/developmental_gradients/data_pmacs/HCPD_xcpd/'
elif dataset == 'HBN':
    datapath = '/cbica/projects/developmental_gradients/data_pmacs/HBN_xcpd/'
elif dataset == 'HCPYA':
    datapath = '/cbica/projects/developmental_gradients/data_pmacs/HCPYA_xcpd/'


# ds = dl.Dataset(datapath)
if dataset == 'HCPYA':
    file_list = glob.glob(os.path.join(datapath, '*.zip'))
else:
    file_list = glob.glob(os.path.join(datapath, 'sub-*.zip'))
file_list.sort()

for f in range(len(file_list)):
    # ds.get(file_list[f])
    # instead of this ^ , I did datalad get in batches
    # (e.g. datalad get sub-0*_xcp-0-3-0.zip)
    # then unzip with this script, then drop

    with zipfile.ZipFile(file_list[f]) as zf:
        fileNames = zf.namelist()

    fileNames = np.array(fileNames)
    if dataset == 'HCPD':
        dataTypes = ['_space-fsLR_atlas-Schaefer417_den-91k_timeseries' +
                     '.ptseries.nii',
                     '_space-fsLR_atlas-Schaefer217_den-91k_timeseries' +
                     '.ptseries.nii',
                     '_space-fsLR_den-91k_qc.csv']
    elif dataset == 'HBN':
        dataTypes = ['_space-fsLR_atlas-Schaefer417_den-91k_timeseries' +
                     '.ptseries.nii',
                     '_space-fsLR_atlas-Schaefer217_den-91k_timeseries' +
                     '.ptseries.nii',
                     '_space-fsLR_den-91k_qc.csv']
    elif dataset == 'HCPYA':
        dataTypes = ['_space-fsLR_seg-4S456Parcels_den-91k_stat' +
                     '-mean_timeseries.ptseries.nii',
                     '_space-fsLR_seg-4S256Parcels_den-91k_stat' +
                     '-mean_timeseries.ptseries.nii',
                     '_motion.tsv']

    wantedFiles = []
    for dataType in dataTypes:
        search = dataType + '+'
        for iFile in fileNames:
            if (re.findall(dataType, iFile)):
                wantedFiles.append(iFile)

    with zipfile.ZipFile(file_list[f], 'r') as zip_ref:
        for specFile in wantedFiles:
            zip_ref.extract(specFile, datapath)
        zip_ref.close()

    # ds.drop(file_list[f])

    print('\nFiles in %s done!' % file_list[f])
