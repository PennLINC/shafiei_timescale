
import os
import glob
import numpy as np
import pandas as pd


dataset = 'HCPYA'
ds_method = 'TRorig'  # can be 'TRmean' or 'TRskip' or 'TRorig' or 'segment'
project_path = '/cbica/projects/developmental_gradients/'

if dataset == 'HCPD':
    inpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
              'timescale/HCPD_acf/%s/' % ds_method)
    fileNames = sorted(glob.glob(inpath +
                                 'sub-*_ses-*_task-*_dir-*_run-*' +
                                 '_space-fsLR_atlas-Schaefer417' +
                                 '_den-91k_timescale_acfsum_%s.npy'
                                 % ds_method))
    inpath_qc = project_path + 'data_pmacs/HCPD_xcpd/xcp_d/'
    fileNames_qc = sorted(glob.glob(inpath_qc +
                                    'sub-*/ses-*/func/sub-*_ses-*_task-*' +
                                    '_dir-*_run-*' +
                                    '_space-fsLR_den-91k_qc.csv'))
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HCPD_acf/%s/' % ds_method)
elif dataset == 'HBN':
    inpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
              'timescale/HBN_acf/%s/' % ds_method)
    fileNames = sorted(glob.glob(inpath +
                                 'sub-*_ses-*_task-*' +
                                 '_space-fsLR_atlas-Schaefer417' +
                                 '_den-91k_timescale_acfsum_%s.npy'
                                 % ds_method))
    inpath_qc = project_path + 'data_pmacs/HBN_xcpd/'
    fileNames_qc = sorted(glob.glob(inpath_qc + 'sub-*/ses-*/func/' +
                                    'sub-*_ses-*_task-*' +
                                    '_space-fsLR_den-91k_qc.csv'))
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HBN_acf/%s/' % ds_method)
elif dataset == 'HCPYA':
    inpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
              'timescale/HCPYA_acf/%s/' % ds_method)
    fileNames = sorted(glob.glob(inpath +
                                 'sub-*_task-*' +
                                 '_space-fsLR_seg-4S456Parcels' +
                                 '_den-91k_stat-mean_timescale_acfsum_%s.npy'
                                 % ds_method))
    inpath_qc = project_path + 'data_pmacs/HCPYA_xcpd/xcpd-0-9-1/'
    fileNames_qc = sorted(glob.glob(inpath_qc + 'sub-*/func/' +
                                    'sub-*_task-*' +
                                    '_dir-*_run-*_motion.tsv'))
    outpath = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'timescale/HCPYA_acf/%s/'
               % ds_method)

subjList = [fileNames[s].split('/')[-1].split('_')[0]
            for s in range(len(fileNames))]
subjList_qc = [fileNames_qc[s].split('/')[-1].split('_')[0]
               for s in range(len(fileNames_qc))]

# check whether subjLists are equal
assert subjList == subjList_qc, 'Subject lists are NOT equal!'

# use the file with the longest filename to generate dataframe columns
check_filename = [len(fileNames[iSubj].split('/')[-1].split('_'))
                  for iSubj in range(len(subjList))]
unique_namelength = np.unique(np.array(check_filename))
maxidx = np.where(np.array(check_filename) == unique_namelength.max())[0][0]

# generate main df
split_name = fileNames[maxidx].split('/')[-1].split('_')
col_names_max = [split_title.split('-')[0] for split_title in split_name[:-3]]
df_main = pd.DataFrame(columns=col_names_max)

# generate main df for qc
if dataset == 'HCPYA':
    subj_qc = pd.read_csv(fileNames_qc[maxidx], delimiter='\t')
else:
    subj_qc = pd.read_csv(fileNames_qc[maxidx])
df_main_qc = pd.DataFrame(columns=list(subj_qc.columns))

# generate empty array for timescale
concat_tau = np.zeros((len(subjList), 400))
df_concat_qc = []
for iSubj in range(len(subjList)):
    subj_tau = np.load(fileNames[iSubj])
    if dataset == 'HCPYA':
        concat_tau[iSubj, :] = subj_tau[:400]
    else:
        concat_tau[iSubj, :] = subj_tau

    split_name = fileNames[iSubj].split('/')[-1].split('_')
    col_names = [split_title.split('-')[0]
                 for split_title in split_name[:-3]]  # -4 for zerocross
    df_temp = pd.DataFrame(columns=col_names)
    col_vals = [split_title.split('-')[1]
                for split_title in split_name[:-3]]  # -4 for zerocross
    df_temp.loc[0] = col_vals

    df_main = pd.concat([df_main, df_temp], ignore_index=True)

    # qc info
    if dataset == 'HCPYA':
        subj_qc = pd.read_csv(fileNames_qc[iSubj], delimiter='\t')
        # Calculate the mean across rows
        mean_series = subj_qc.mean(axis=0)
        # Convert the mean Series to a dataframe with one row + reset index
        subj_qc_mean = pd.DataFrame(mean_series).T
        subj_qc_mean = subj_qc_mean.reset_index(drop=True)

        df_main_qc = pd.concat([df_main_qc, subj_qc_mean], ignore_index=True)
    else:
        subj_qc = pd.read_csv(fileNames_qc[iSubj])
        df_main_qc = pd.concat([df_main_qc, subj_qc], ignore_index=True)


# change into array: in 17-network ordering (all data except for HCPYA)
concat_tau = np.array(concat_tau)

reorder_dir = (project_path + 'gitrepo/shafiei_timescale/data/' +
               'SchaeferParcellation/schaefer_ordering_mapper/')
mapped = pd.read_csv(reorder_dir + 'Schaefer_400-17_mappedto_400-7.csv')

if dataset == 'HCPYA':
    # HCPYA is already in 7-network ordering
    # create dataframe
    final_df = pd.DataFrame(concat_tau, columns=list(mapped['output_roi']))
    final_df = pd.concat([final_df, df_main], axis=1)
else:
    # matrix permuted to 7-network ordering
    concat_tau_17to7 = concat_tau[:, mapped['mapped_indices'].values]

    # create dataframe
    final_df = pd.DataFrame(concat_tau_17to7,
                            columns=list(mapped['output_roi']))
    final_df = pd.concat([final_df, df_main], axis=1)

# save files
if not os.path.exists(os.path.join(outpath, 'concat/')):
    os.makedirs(os.path.join(outpath, 'concat/'))

final_df.to_csv(outpath + 'concat/' +
                '%s_concat_timescale_acf_%s_Schaefer_400-7.csv'
                % (dataset, ds_method), index=False)
df_main_qc.to_csv(outpath + 'concat/%s_concat_functional_qc.csv' % dataset,
                  index=False)
