import os
import numpy as np
import pandas as pd

###################
# HCPYA dataset
###################
# # Sample Selection process
# 1) include passing T1 QC: all scans in dataset have survived T1 QC already
# 2) include meanFD < 0.2


dataset = 'HCPYA'
ds_method = 'TRorig'
task = 'rest'
dtype = 'timescale'  # 'dtype' is timescale for this project
if dtype == 'timescale':
    metric = 'timescale_acf'
else:
    metric = dtype

datapath = ('/cbica/projects/developmental_gradients/' +
            'gitrepo/shafiei_timescale/data/')

# load timescale data
ts_data = pd.read_csv(datapath + '%s/%s_acf/%s/concat/'
                      % (dtype, dataset, ds_method) +
                      '%s_concat_%s_%s_Schaefer_400-7.csv'
                      % (dataset, metric, ds_method),
                      low_memory=False)

###################
# load demographics
###################
# load demographics data
df_demogs = pd.read_csv(datapath + 'qc_demos/demo/'
                        'RESTRICTED_gshafiei_12_1_2022_10_58_3.csv')
df_demogs_unrestric = pd.read_csv(datapath + 'qc_demos/demo/unrestricted_' +
                                  'gshafiei_9_18_2024_11_42_52.csv')
df_demogs['sex'] = df_demogs_unrestric['Gender'].copy()
###################
# fc qc
###################
# load fc qc file: they are with timescale data
fc_qc = pd.read_csv(datapath + '%s/%s_acf/%s/concat/'
                    % (dtype, dataset, ds_method) +
                    '%s_concat_functional_qc.csv' % (dataset),
                    low_memory=False)

# add task and sub id info
fc_qc['task'] = ts_data['task'].copy()
fc_qc['sub'] = ts_data['sub'].copy()

# get task specific data
fc_qc_task = fc_qc[fc_qc['task'].str.contains('%s' % task)]
fc_qc_task.reset_index(inplace=True, drop=True)
# set datatype for subj ids in fc qc file
fc_qc_task['sub'] = fc_qc_task['sub'].astype(str)

###################
# timescale data + fc qc + add demogs
###################
# get task specific data
ts_data_task = ts_data[ts_data['task'].str.contains('%s' % task)]
ts_data_task.reset_index(inplace=True, drop=True)
# set datatype for subj ids in ts_data file
ts_data_task['sub'] = ts_data_task['sub'].astype(str)

# add fc qc info to ts data
if list(ts_data_task['sub']) == list(fc_qc_task['sub']):
    ts_data_task['meanFD'] = fc_qc_task['framewise_displacement']
else:
    raise ValueError('fc qc and ts subj lists are NOT equal!')

print('\nInitial sample size for HCYPA rest: %s'
      % ts_data_task['sub'].unique().shape)

# remove meanFD < 0.2
ts_data_task = ts_data_task.query('meanFD < 0.2')
ts_data_task.reset_index(inplace=True, drop=True)

# Only average the numeric columns to avoid TypeError
numeric_data = ts_data_task.select_dtypes(include='number')
avg_runs = numeric_data.groupby(ts_data_task['sub']).mean()

# Keep only the first 400 columns and 'meanFD' (if it exists)
col_keep = list(avg_runs.columns[:400])
if 'meanFD' in avg_runs.columns and 'meanFD' not in col_keep:
    col_keep.append('meanFD')

avg_runs = avg_runs[col_keep]

# Add 'sub-' prefix to the subject IDs
ts_sub_list = list(avg_runs.index.values)
updated_subj_list = ['sub-' + str(iSubj) for iSubj in ts_sub_list]
avg_runs['participant_id'] = updated_subj_list
avg_runs.reset_index(inplace=True, drop=True)

print('\nSample size after motion exclusion: %s'
      % avg_runs['participant_id'].unique().shape)

# Add 'sub-' prefix to the subject IDs in demogs info as well
demogs_sub_list = list(df_demogs['Subject'].values)
updated_subj_list = ['sub-' + str(iSubj) for iSubj in demogs_sub_list]
df_demogs['participant_id'] = updated_subj_list

# find intersect of ids in avg runs and demogs
x = np.array(avg_runs['participant_id'])
y = np.array(df_demogs['participant_id'])
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

avg_runs_wdemos = avg_runs.iloc[x_ind, :]
avg_runs_wdemos.reset_index(inplace=True, drop=True)

df_demogs_shared = df_demogs.iloc[y_ind, :]
df_demogs_shared.reset_index(inplace=True, drop=True)

# add demographics info
final_df = pd.merge(avg_runs_wdemos, df_demogs_shared,
                    on='participant_id')

# add mean TS and std TS
meanTS = np.nanmean(final_df.iloc[:, 0:400], axis=1)
stdTS = np.nanstd(final_df.iloc[:, 0:400], axis=1)

if dtype == 'timescale':
    mean_colname = 'meanTS'
    std_colname = 'stdTS'

final_df.insert(loc=400, column=mean_colname, value=meanTS)
final_df.insert(loc=401, column=std_colname, value=stdTS)

# age in years is from restricted demogs
final_df['age'] = final_df['Age_in_Yrs'].copy()
final_df['study'] = ['HCPYA']*len(final_df)

outpath = datapath + '%s/%s_%s_%s_concat/' % (dtype, dataset,
                                              task, ds_method)
if not os.path.exists(outpath):
    os.makedirs(outpath)

# save output
final_df.to_csv(outpath +
                '%s_%s_concat_%s_%s_Schaefer_400-7_forR.tsv'
                % (dataset, task, metric, ds_method),
                index=False, sep='\t')
