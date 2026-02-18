import os
import time
import glob
import numpy as np
import scipy.stats
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


#################################
def autocorr_decay(dk, A, tau, B):
    return A*(np.exp(-(dk/tau))+B)
#################################


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

        TR = TR_orig
        ts_data_zscore = scipy.stats.zscore(ts_data, axis=0).T

        nParcels = np.shape(ts_data_zscore)[0]
        maxlag = 100
        taus = np.zeros((nParcels,))
        As = np.zeros((nParcels,))
        Bs = np.zeros((nParcels,))
        xdata = np.arange(maxlag)

        for roi in range(nParcels):
            tscale = np.correlate(ts_data_zscore[roi, :],
                                  ts_data_zscore[roi, :],
                                  mode='full')
            start = len(tscale)//2 + 1
            # Normalize
            tscale = np.divide(tscale, np.max(tscale))
            stop = start + maxlag
            tscale = tscale[start:stop]
            # Sometimes curve_fit cannot find a good solution, so if it can't,
            # enter a 'nan' for this subject's ROI
            try:
                A, tau, B = scipy.optimize.curve_fit(autocorr_decay, xdata,
                                                     tscale,
                                                     p0=[0,
                                                         np.random.rand(1)[0] +
                                                         0.01, 0],
                                                     bounds=(([0, 0, -np.inf],
                                                              [np.inf, np.inf,
                                                               np.inf])),
                                                     method='trf')[0]
                As[roi,] = A
                taus[roi,] = tau
                Bs[roi,] = B
            except:
                print('ROI', roi, 'NaN')
                As[roi,] = np.nan
                taus[roi,] = np.nan
                Bs[roi,] = np.nan

        new_fileName = dataFile.split('/')[-1].replace('timeseries.ptseries' +
                                                       '.nii',
                                                       'timescale_tau_' +
                                                       '%s.npy' % ds_method)
        np.save(outpath + new_fileName, taus)

    tend = time.time()
    print('\nSubj', iFile, 'of', len(fileNames), 'done!',
          '\nRunning time = ', tend-tstart, 'seconds!')
