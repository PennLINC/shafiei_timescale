import os
import numpy as np
import pandas as pd

###################
# HCPD dataset
###################
# # Sample Selection process
# 1) exclude participants with medical conditions affecting brain function,
# gross neurological abnormalities
# 2) include passing T1 QC: all scans in dataset have survived T1 QC already
# 3) include meanFD < 0.2
# 4) include individuals within the age range of 8 to 22 years old


dataset = 'HCPD'
ds_method = 'TRorig'
task = 'rest'
dtype = 'timescale'
if dtype == 'timescale':  # 'dtype' is timescale for this project
    metric = 'timescale_acf'
else:
    metric = dtype

datapath = ('/cbica/projects/developmental_gradients/' +
            'gitrepo/shafiei_timescale/data/')

# load timescale data
ts_data = pd.read_csv(datapath + '%s/%s_acf/%s/concat/'
                      % (dtype, dataset, ds_method) +
                      '%s_concat_%s_%s_Schaefer_200-7.csv'
                      % (dataset, metric, ds_method),
                      low_memory=False)

###################
# health exclusion
###################
# load hcpd_demographics file that includes health data
health_data = pd.read_csv(datapath + 'qc_demos/demo/hcpd_demographics.csv',
                          low_memory=False)

# cancer/leukemia: HCD0641341
cancer_idx = np.where(health_data['medhis_2e'] == 1)[0]
# lead poisoning: HCD0392144, HCD0641341
lead_idx = np.where(health_data['medhis_2k'] == 1)[0]
# sickle cell anemia: HCD0641341
sickle_idx = np.where(health_data['medhis_2p'] == 1)[0]
# accidental poisoning: HCD0529751, HCD2855875
poisoning_idx = np.where(health_data['medhis_6q'] == 1)[0]
# multiple sclerosis: HCD0641341, HCD1617450
ms_idx = np.where(health_data['ms'] == 1)[0]
# seizure: ph9: HCD0360838, HCD0478356; cfmh: HCD0641341, HCD1277149
seizure_ph9_idx = np.where(health_data['ph_9'] == 1)[0]
epilepsy_cfmh_idx = np.where(health_data['cfmh_chd_seizure'] == 1)[0]
# brain injury: 14 participants: HCD0041822, HCD0514738, HCD0641341,
# HCD1162334, HCD1534244, HCD1785572, HCD1801140, HCD1886477, HCD1903350,
# HCD2181848, HCD2256651, HCD2378463, HCD2496570, HCD2741961
brain_injury_idx = np.where(health_data['seq1c_2'] == 1)[0]

# remove them all from health data
all_medexclude = (list(cancer_idx) + list(lead_idx) + list(sickle_idx) +
                  list(poisoning_idx) + list(ms_idx) + list(seizure_ph9_idx) +
                  list(epilepsy_cfmh_idx) + list(brain_injury_idx))
health_data.drop(np.unique(np.array(all_medexclude)), inplace=True)
health_data.reset_index(inplace=True, drop=True)

###################
# load demographics and exclude health excluded participants
###################
# load demographics data
# note: this file comes from
# https://github.com/PennLINC/RBC_demo_pheno/tree/main/UPDATED_TSVS
df_demogs = pd.read_csv(datapath + 'qc_demos/demo/'
                        '%s_participants.tsv' % dataset.lower(),
                        delimiter='\t')

# replace 'HCD' with 'sub-' in health data subject id
health_subj_list = list(health_data['src_subject_id'])
health_subj_list = [subjid.replace('HCD', 'sub-')
                    for subjid in health_subj_list]

# find intersect of health and demogs ids and regenerate the lists
x = np.array(df_demogs['participant_id'])
y = np.array(health_subj_list)
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

df_demogs_health = df_demogs.iloc[x_ind, :]
df_demogs_health.reset_index(inplace=True, drop=True)

###################
# fc qc
###################
# load fc qc file: they are with timescale data
fc_qc = pd.read_csv(datapath + '%s/%s_acf/%s/concat/'
                    % (dtype, dataset, ds_method) +
                    '%s_concat_functional_qc.csv' % (dataset),
                    low_memory=False)
# get task specific data
fc_qc_task = fc_qc[fc_qc['task'].str.contains('%s' % task)]
fc_qc_task.reset_index(inplace=True, drop=True)
# zeropad subj ids in fc qc file
fc_qc_task['sub'] = fc_qc_task['sub'].astype(str)
fc_qc_task['sub'] = fc_qc_task['sub'].str.zfill(7)

###################
# timescale data + fc qc + add demogs
###################
# get task specific data
ts_data_task = ts_data[ts_data['task'].str.contains('%s' % task)]
ts_data_task.reset_index(inplace=True, drop=True)
# zeropad subj ids in ts_data file
ts_data_task['sub'] = ts_data_task['sub'].astype(str)
ts_data_task['sub'] = ts_data_task['sub'].str.zfill(7)

# add fc qc info to ts data
if list(ts_data_task['sub']) == list(fc_qc_task['sub']):
    ts_data_task['meanFD'] = fc_qc_task['meanFD']
    ts_data_task['nVolsRemoved'] = fc_qc_task['nVolsRemoved']
    ts_data_task['num_censored_volumes'] = fc_qc_task['num_censored_volumes']
else:
    raise ValueError('fc qc and ts subj lists are NOT equal!')

print('\nInitial sample size for HCPD rest: %s'
      % ts_data_task['sub'].unique().shape)

# remove meanFD < 0.2
ts_data_task = ts_data_task.query('meanFD < 0.2')
ts_data_task.reset_index(inplace=True, drop=True)

# check whether any volumes were removed or censored: drop them
check_vol_removed = np.where(ts_data_task['nVolsRemoved'] >= 1)[0].any()
check_vol_censored = np.where(ts_data_task['num_censored_volumes']
                              >= 1)[0].any()
if check_vol_removed or check_vol_censored:
    vol_removed_idx = np.where(ts_data_task['nVolsRemoved'] >= 1)[0]
    vol_censored_idx = np.where(ts_data_task['num_censored_volumes'] >= 1)[0]
    all_vol_idx = list(vol_removed_idx) + list(vol_censored_idx)
    ts_data_task.drop(np.array(all_vol_idx), inplace=True)
    ts_data_task.reset_index(inplace=True, drop=True)

# Only average the numeric columns to avoid TypeError
numeric_data = ts_data_task.select_dtypes(include='number')
avg_runs = numeric_data.groupby(ts_data_task['sub']).mean()

# Keep only the first 200 columns and 'meanFD' (if it exists)
col_keep = list(avg_runs.columns[:200])
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

# find intersect of ids in avg runs and demogs
x = np.array(avg_runs['participant_id'])
y = np.array(df_demogs_health['participant_id'])
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

avg_runs_wdemos = avg_runs.iloc[x_ind, :]
avg_runs_wdemos.reset_index(inplace=True, drop=True)

df_demogs_shared = df_demogs_health.iloc[y_ind, :]
df_demogs_shared.reset_index(inplace=True, drop=True)

# add demographics info to TS dataframe (automatically excludes
# medical exclusions)
# then add mean TS and std TS + demogs
final_df = pd.merge(avg_runs_wdemos, df_demogs_shared,
                    on='participant_id')
meanTS = np.nanmean(final_df.iloc[:, 0:200], axis=1)
stdTS = np.nanstd(final_df.iloc[:, 0:200], axis=1)

if dtype == 'timescale':
    mean_colname = 'meanTS'
    std_colname = 'stdTS'

final_df.insert(loc=200, column=mean_colname, value=meanTS)
final_df.insert(loc=201, column=std_colname, value=stdTS)

print('\nSample size after medical exclusion: %s'
      % final_df['participant_id'].unique().shape)

outpath = datapath + '%s/%s_%s_%s_concat/' % (dtype, dataset,
                                              task, ds_method)
if not os.path.exists(outpath):
    os.makedirs(outpath)

# save output
final_df.to_csv(outpath +
                '%s_%s_concat_%s_%s_Schaefer_200-7_forR.tsv'
                % (dataset, task, metric, ds_method),
                index=False, sep='\t')

###################
# save 8-22 age range
###################
temp = final_df.copy()
final_df_ageabove8 = temp.query('age >= 8')
final_df_ageabove8.reset_index(inplace=True, drop=True)

final_df_age822 = final_df_ageabove8.query('age <= 22')
final_df_age822.reset_index(inplace=True, drop=True)

print('\nFinal sample size for HCPD rest: %s'
      % final_df_age822['participant_id'].unique().shape)

outpath = datapath + '%s/%s_%s_age_%s_concat/' % (dtype, dataset,
                                                  task, ds_method)
if not os.path.exists(outpath):
    os.makedirs(outpath)

final_df_age822.to_csv(outpath + '%s_%s_age_concat_%s_%s_Schaefer_200-7'
                       % (dataset, task, metric, ds_method) +
                       '_forR.tsv',
                       index=False, sep='\t')
