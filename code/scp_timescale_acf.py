import os
import time
import glob
import random
import numpy as np
import scipy.stats
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf

# set paths
dataset = 'HCPYA'
ds_method = 'TRorig'  # can be 'TRmean' or 'TRskip' or 'TRorig' or 'segment'
project_path = '/cbica/projects/developmental_gradients/'

if dataset == 'HCPD':
    inpath = project_path + 'data_pmacs/HCPD_xcpd/xcp_d/'
    fileNames = glob.glob(inpath + 'sub-*/ses-*/func/' +
                          'sub-*_ses-*_task-*_dir-*_run-*_space-fsLR_' +
                          'atlas-Schaefer417_den-91k_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HCPD_acf/%s/' % ds_method)
    TR_orig = 0.8  # s
elif dataset == 'HBN':
    inpath = project_path + 'data_pmacs/HBN_xcpd/'
    fileNames = glob.glob(inpath + 'sub-*/ses-*/func/' +
                          'sub-*_ses-*_task-*_space-fsLR_' +
                          'atlas-Schaefer417_den-91k_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HBN_acf/%s/' % ds_method)
elif dataset == 'HCPYA':
    inpath = project_path + 'data_pmacs/HCPYA_xcpd/xcpd-0-9-1/'
    fileNames = glob.glob(inpath + 'sub-*/func/' +
                          'sub-*_task-*_dir-*_run-*_space-fsLR_' +
                          'seg-4S456Parcels_' +
                          'den-91k_stat-mean_timeseries.ptseries.nii')
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               '/timescale/HCPYA_acf/%s/' % ds_method)
    TR_orig = 0.72  # s

# make outpath directory if it doesn't exist
if not os.path.exists(outpath):
    os.makedirs(outpath)

# to randomly plot acf of some individuals
random.seed(32)
randidx = random.sample(range(len(fileNames)), 50)

# to plot specific rois
reorder_dir = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'SchaeferParcellation/schaefer_ordering_mapper/')
mapped = pd.read_csv(reorder_dir + 'Schaefer_400-17_mappedto_400-7.csv')

# sa axis rank for 7 network order
sa_rank_7 = pd.read_csv(project_path + 'gitrepo/shafiei_timescale/data/' +
                        'SchaeferParcellation/' +
                        'SArank_schaefer400_7Networks.csv')

# want to plot 8 acfs for regions with sa rank of 15, 30, 45, 60, 75, 90
sa_region = list(np.linspace(1, 400, 15, dtype=int))

# list of regions to plot in 7rsn ordering
sa_region_idx_in7 = [np.where(sa_rank_7['SArank'].values == region)[0][0]
                     for region in sa_region]

# find these regions in 17rsn ordering; so where would each region be
# in 17 rsn ordering mapping
sa_region_idx_in17 = mapped['mapped_indices'].iloc[sa_region_idx_in7]


for iFile, dataFile in enumerate(fileNames):
    tstart = time.time()

    outFile = outpath + dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                        '.nii',
                                                        'timescale_acfsum_' +
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
            ts_data_zscore = scipy.stats.zscore(ts_downsampled, axis=0).T
            plot_nlags = 40
        elif ds_method == 'TRskip':
            # or select every 4th point
            ts_copy = ts_data.copy()
            ts_downsampled = ts_copy[0::4, :]
            TR = TR_orig * 4
            ts_data_zscore = scipy.stats.zscore(ts_downsampled, axis=0).T
            plot_nlags = 40
        elif ds_method == 'TRorig':
            # orig TR
            TR = TR_orig
            ts_data_zscore = scipy.stats.zscore(ts_data, axis=0).T
            plot_nlags = 50
        elif ds_method == 'segment':
            # orig TR, short segment
            TR = TR_orig
            ts_data_segment = ts_data[:150, :]
            ts_data_zscore = scipy.stats.zscore(ts_data_segment, axis=0).T
            plot_nlags = 40

        nParcels = np.shape(ts_data_zscore)[0]
        nTimepoints = np.shape(ts_data_zscore)[1]
        acf_sum = np.zeros((nParcels,))

        for roi in range(nParcels):
            acf_all = acf(ts_data_zscore[roi, :], nlags=nTimepoints-1)

            try:
                negidx = np.where(acf_all < 0)[0]
                first_zerocross = np.min(negidx)
                acf_sum[roi,] = np.sum(acf_all[:first_zerocross]) * TR
            except Exception:
                print('ROI', roi, 'infinite')
                acf_sum[roi,] = np.nan

        new_fileName = dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                       '.nii',
                                                       'timescale_acfsum_' +
                                                       '%s.npy' % ds_method)
        np.save(outpath + new_fileName, acf_sum)

        if iFile in randidx:
            if dataset != 'HCPYA':
                toplot_regions = sa_region_idx_in17.values
            else:
                toplot_regions = sa_region_idx_in7
            plt.ion()
            fig, ax = plt.subplots(3, 5, figsize=(15, 9))
            for r, roi in enumerate(toplot_regions):
                if r < 5:
                    plot_acf(ts_data_zscore[roi, :], lags=plot_nlags,
                             ax=ax[0][r],
                             title='SArank: %s' % str(sa_region[r]))
                elif r >= 5 and r < 10:
                    plot_acf(ts_data_zscore[roi, :], lags=plot_nlags,
                             ax=ax[1][r-5],
                             title='SArank: %s' % str(sa_region[r]))
                elif r >= 10:
                    plot_acf(ts_data_zscore[roi, :], lags=plot_nlags,
                             ax=ax[2][r-10],
                             title='SArank: %s' % str(sa_region[r]))
            figName = dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                      '.nii',
                                                      'timescale_acf_' +
                                                      '%s.png' % ds_method)
            if not os.path.exists(os.path.join(outpath, 'png/')):
                os.makedirs(os.path.join(outpath, 'png/'))
            plt.savefig(outpath + 'png/' + figName)
            plt.close()
    tend = time.time()
    print('\nSubj', iFile, 'of', len(fileNames), 'done!',
          '\nRunning time = ', tend-tstart, 'seconds!')
