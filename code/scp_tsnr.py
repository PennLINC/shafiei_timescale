import os
import time
import glob
import numpy as np
import pandas as pd
import nibabel as nib

# set paths
dataset = 'HCPD'
ds_method = 'TRorig'  # can be 'TRmean' or 'TRskip' or 'TRorig' or 'segment'
project_path = '/cbica/projects/developmental_gradients/'

if dataset == 'HCPD':
    inpath = project_path + 'data_pmacs/HCPD_xcpd/xcp_d/'
    fileNames = glob.glob(inpath + 'sub-*/ses-*/func/' +
                          'sub-*_ses-*_task-*_dir-*_run-*_space-fsLR_' +
                          'atlas-Schaefer417_den-91k_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HCPD_tsnr/%s/' % ds_method)
    TR_orig = 0.8  # s
elif dataset == 'HBN':
    inpath = project_path + 'data_pmacs/HBN_xcpd/'
    fileNames = glob.glob(inpath + 'sub-*/ses-*/func/' +
                          'sub-*_ses-*_task-*_space-fsLR_' +
                          'atlas-Schaefer417_den-91k_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HBN_tsnr/%s/' % ds_method)
elif dataset == 'HCPYA':
    inpath = project_path + 'data_pmacs/HCPYA_xcpd/xcpd-0-9-1/'
    fileNames = glob.glob(inpath + 'sub-*/func/' +
                          'sub-*_task-*_dir-*_run-*_space-fsLR_' +
                          'seg-4S456Parcels_' +
                          'den-91k_stat-mean_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               '/timescale/HCPYA_tsnr/%s/' % ds_method)
    TR_orig = 0.72  # s

# make outpath directory if it doesn't exist
if not os.path.exists(outpath):
    os.makedirs(outpath)

for iFile, dataFile in enumerate(fileNames):
    tstart = time.time()

    outFile = outpath + dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                        '.nii',
                                                        'tsnr_' +
                                                        '%s.npy' % ds_method)

    if dataset == 'HBN':
        ses_info = dataFile.split('/')[7]
        if ses_info == 'ses-HBNsiteSI':
            TR_orig = 1.45  # s
        else:
            TR_orig = 0.8  # s

    if not os.path.isfile(outFile):
        ts_file = nib.load(dataFile)
        ts_data = np.array(ts_file.get_fdata())

        # orig TR or short TR: want to downsample from 0.8s TR to 3s TR:
        if ds_method == 'TRmean':
            # so 3/0.8 = 3.75 --> 4 (avergae every 4 points)
            ts_df = pd.DataFrame(ts_data)
            ts_downsampled_df = ts_df.groupby(np.arange(len(ts_df))//4).mean()
            ts_downsampled = np.array(ts_downsampled_df)
            TR = TR_orig * 4
            ts_data_transpose = ts_downsampled.T
            plot_nlags = 40
        elif ds_method == 'TRskip':
            # or select every 4th point
            ts_copy = ts_data.copy()
            ts_downsampled = ts_copy[0::4, :]
            TR = TR_orig * 4
            ts_data_transpose = ts_downsampled.T
            plot_nlags = 40
        elif ds_method == 'TRorig':
            # orig TR
            TR = TR_orig
            ts_data_transpose = ts_data.T
            plot_nlags = 50
        elif ds_method == 'segment':
            # orig TR, short segment
            TR = TR_orig
            ts_data_segment = ts_data[:150, :]
            ts_data_transpose = ts_data_segment.T
            plot_nlags = 40

        nParcels = np.shape(ts_data_transpose)[0]
        nTimepoints = np.shape(ts_data_transpose)[1]
        tsnr = np.zeros((nParcels,))

        for roi in range(nParcels):
            tsnr[roi,] = np.mean(ts_data_transpose[roi, :])/np.std(
                ts_data_transpose[roi, :])

        new_fileName = dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                       '.nii',
                                                       'tsnr_' +
                                                       '%s.npy' % ds_method)
        np.save(outpath + new_fileName, tsnr)

    tend = time.time()
    print('\nSubj', iFile, 'of', len(fileNames), 'done!',
          '\nRunning time = ', tend-tstart, 'seconds!')
