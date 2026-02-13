import os
import numpy as np
import pandas as pd

###################
# HBN dataset
###################
# # Sample Selection process
# 1) include passing T1 QC
# 2) include meanFD < 0.2
# 3) include indiviudals within the age range of 8 to 22 years old

dataset = 'HBN'
ds_method = 'TRorig'
task = 'rest'
filter_qc = True
dtype = 'timescale'
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

# update how subject IDs (to be compatible with demogs later)
# add sub- to the beginning of ids
ts_sub_list = list(ts_data['sub'].values)
updated_subj_list = ['sub-' + str(iSubj) for iSubj in ts_sub_list]
ts_data['participant_id'] = updated_subj_list
ts_data.drop(columns=['sub'], inplace=True)

###################
# load demographics
###################
# load demographics data
# note: this file comes from
# https://github.com/PennLINC/RBC_demo_pheno/tree/main/+
# UPDATED_TSV_WITH_ETHNICITY_BINNED.zip
df_demogs = pd.read_csv(datapath + 'qc_demos/demo/' +
                        '%s_participants.tsv' % dataset.lower(),
                        delimiter='\t')

# add study site from demogs to ts_data
ts_data = ts_data.merge(df_demogs[['participant_id', 'study_site']],
                        on='participant_id',
                        how='left')

###################
# drop SI site from timescale data + keep task-specific data
###################
ts_data = ts_data[~ts_data['study_site'].str.contains('HBNsiteSI')]
ts_data.reset_index(inplace=True, drop=True)

# keep task specific data
ts_data_task = ts_data[ts_data['task'].str.contains('%s' % task)]
ts_data_task.reset_index(inplace=True, drop=True)

print('\nInitial sample size for HBN rest: %s'
      % ts_data_task['participant_id'].unique().shape)

###################
# T1 qc
###################
# load T1 qc
# note this file comes from Freesurfer data on RBC:
# https://github.com/ReproBrainChart/HBN_FreeSurfer
t1qc = pd.read_csv(datapath + 'qc_demos/T1QC/study-%s_desc-T1_qc.tsv'
                   % dataset, low_memory=False, delimiter='\t')

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

# update how subject IDs (to be compatible with demogs later)
# add sub- to the beginning of ids
fcqc_sub_list = list(fc_qc_task['sub'].values)
updated_subj_list = ['sub-' + str(iSubj) for iSubj in fcqc_sub_list]
fc_qc_task['participant_id'] = updated_subj_list
fc_qc_task.drop(columns=['sub'], inplace=True)

# drop SI from fc qc data
fc_qc_task = fc_qc_task[~fc_qc_task['ses'].str.contains('HBNsiteSI')]
fc_qc_task.reset_index(inplace=True, drop=True)

# add fc qc info to ts data
if (list(ts_data_task['participant_id']) == list(fc_qc_task['participant_id'])
        and list(ts_data_task['study_site']) == list(fc_qc_task['ses'])):
    ts_data_task['meanFD'] = fc_qc_task['meanFD']
    ts_data_task['nVolsRemoved'] = fc_qc_task['nVolsRemoved']
    ts_data_task['num_censored_volumes'] = fc_qc_task['num_censored_volumes']
else:
    raise ValueError('fc qc and ts subj lists are NOT equal!')

###################
# timescale data + T1qc
###################
# add T1 QC info (euler and qc determination label) from t1qc to ts_data_task
ts_data_task = ts_data_task.merge(t1qc[['participant_id',
                                        'euler', 'qc_determination']],
                                  on='participant_id', how='left')

# exclude individuals with failed T1
if filter_qc:
    ts_data_task = ts_data_task[
        ~ts_data_task['qc_determination'].str.contains('Fail')]
    ts_data_task.reset_index(inplace=True, drop=True)

print('\nSample size after T1 QC exclusion: %s'
      % ts_data_task['participant_id'].unique().shape)

###################
# timescale data: apply fc qc + add demogs
###################
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

print('\nSample size after motion exclusion: %s'
      % ts_data_task['participant_id'].unique().shape)

# groupby sub only
avg_runs = ts_data_task.groupby('participant_id').mean(numeric_only=True)
col_keep = list(avg_runs.columns[:400])
col_keep.extend(['meanFD'])
avg_runs = avg_runs[col_keep]
avg_runs['participant_id'] = list(avg_runs.index.values)
avg_runs.reset_index(inplace=True, drop=True)

# add T1qc info back to avg_runs
avg_runs = avg_runs.merge(t1qc[['participant_id',
                                'euler', 'qc_determination']],
                          on='participant_id', how='left')

# add demogs to avg_runs
avg_runs = avg_runs.merge(df_demogs, on='participant_id', how='left')

# add mean TS and std TS
final_df = avg_runs.copy()
meanTS = np.nanmean(final_df.iloc[:, 0:400], axis=1)
stdTS = np.nanstd(final_df.iloc[:, 0:400], axis=1)

if dtype == 'timescale':
    mean_colname = 'meanTS'
    std_colname = 'stdTS'

final_df.insert(loc=400, column=mean_colname, value=meanTS)
final_df.insert(loc=401, column=std_colname, value=stdTS)

# save output
outpath = datapath + '%s/%s_%s_noSI_%s_concat/' % (dtype, dataset,
                                                   task, ds_method)

if not os.path.exists(outpath):
    os.makedirs(outpath)

final_df.to_csv(outpath + '%s_%s_noSI_concat_%s_%s_Schaefer_400-7_forR.tsv'
                % (dataset, task, metric, ds_method),
                index=False, sep='\t')

###################
# save for 8-22 age range
###################
# age 8-22
final_df_ageabove8 = final_df.query('age >= 8')
final_df_ageabove8.reset_index(inplace=True, drop=True)

final_df_age822 = final_df_ageabove8.query('age <= 22')
final_df_age822.reset_index(inplace=True, drop=True)

print('\nSample size after motion exclusion: %s'
      % final_df_age822['participant_id'].unique().shape)

outpath = datapath + '%s/%s_%s_noSI_age_%s_concat/' % (dtype, dataset,
                                                       task, ds_method)
if not os.path.exists(outpath):
    os.makedirs(outpath)

final_df_age822.to_csv(outpath + '%s_%s_noSI_age_concat_%s_%s_Schaefer_400-7'
                       % (dataset, task, metric, ds_method) + '_forR.tsv',
                       index=False, sep='\t')

###################
# save 600 low motion
###################
final_df_age822_sorted = final_df_age822.sort_values('meanFD')
final_df_age822_sorted.reset_index(inplace=True, drop=True)
final_df_age822_sorted = final_df_age822_sorted.iloc[:600, :]

outpath = datapath + '%s/%s_%s_noSI_age_600LM_%s_concat/' % (dtype, dataset,
                                                             task, ds_method)

if not os.path.exists(outpath):
    os.makedirs(outpath)

final_df_age822_sorted.to_csv((outpath + ('%s_%s_noSI_age_600LM_concat_%s_%s' +
                               '_Schaefer_400-7')
                               % (dataset, task, metric, ds_method) +
                               '_forR.tsv'), index=False, sep='\t')

###################
# save 600 low pfactor
###################
final_df_age822_sorted = final_df_age822.sort_values('p_factor_mcelroy' +
                                                     '_harmonized_all_samples')
final_df_age822_sorted.reset_index(inplace=True, drop=True)
final_df_age822_sorted = final_df_age822_sorted.iloc[:600, :]

outpath = datapath + '%s/%s_%s_noSI_age_600LP_%s_concat/' % (dtype, dataset,
                                                             task, ds_method)

if not os.path.exists(outpath):
    os.makedirs(outpath)

final_df_age822_sorted.to_csv((outpath + ('%s_%s_noSI_age_600LP_concat_%s_%s' +
                               '_Schaefer_400-7')
                               % (dataset, task, metric, ds_method) +
                               '_forR.tsv'), index=False, sep='\t')
