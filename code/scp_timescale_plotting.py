################
# plotting script for timescale
################
import scipy
import numpy as np
import pandas as pd
import fcn_timescale
import seaborn as sns
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# import matplotlib
# matplotlib.use('Agg')  # non-interactive backend
# import matplotlib.pyplot as plt

# set paths
inpath = '/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/'
parcellationDir = inpath + 'data/SchaeferParcellation/'

# atlas data
lhlabels = (parcellationDir + 'fslr32k/' +
            'Schaefer2018_400Parcels_7Networks_order_lh.label.gii')
rhlabels = (parcellationDir + 'fslr32k/' +
            'Schaefer2018_400Parcels_7Networks_order_rh.label.gii')
labelinfo = np.loadtxt(parcellationDir + 'fslr32k/' +
                       'Schaefer2018_400Parcels_7Networks_order_info.txt',
                       dtype='str', delimiter='\t')
regionlabels = []
for row in range(0, len(labelinfo), 2):
    regionlabels.append(labelinfo[row])

# surface data
surf_path = inpath + 'data/surfaces/'

surfaces = [surf_path + 'L.sphere.32k_fs_LR.surf.gii',
            surf_path + 'R.sphere.32k_fs_LR.surf.gii']

# load sa rank
sa_rank = pd.read_csv(parcellationDir + 'SArank_schaefer400_7Networks.csv')

# get custom colormaps
# cmap_seq, cmap_seq_r, megcmap, megcmap_r, megcmap2, categ_cmap
_, _, megcmap, _, _, _, _, _, cmap_YlPu = fcn_timescale.make_colormaps()

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Myriad Pro']
plt.rcParams['font.size'] = 18.0

################
# correlate maps themselves
################
hbn_ts = pd.read_csv(inpath + 'data/timescale/' +
                     'HBN_rest_noSI_age_TRorig_concat/' +
                     'HBN_rest_noSI_age_concat_timescale_acf_TRorig' +
                     '_Schaefer_400-7_forR.tsv', sep='\t')
hcpd_ts = pd.read_csv(inpath + 'data/timescale/HCPD_rest_age_TRorig_concat/' +
                      'HCPD_rest_age_concat_timescale_acf_TRorig' +
                      '_Schaefer_400-7_forR.tsv', sep='\t')

x = pd.DataFrame(np.nanmean(np.array(hcpd_ts.iloc[:, :400]), axis=0))
y = pd.DataFrame(np.nanmean(np.array(hbn_ts.iloc[:, :400]), axis=0))

combined_df = pd.concat([x, y], axis=1)
labels = ['hcpd', 'hbn']
combined_df.columns = labels

# spin p
nspins = 10000

x = combined_df['hcpd'].values
y = combined_df['hbn'].values

corrval = scipy.stats.spearmanr(x, y)[0]
pval = fcn_timescale.get_spinp(x, y, corrval=corrval, nspin=nspins,
                               lhannot=lhlabels, rhannot=rhlabels,
                               corrtype='spearman',
                               surfpath=surf_path)

plt.ion()
title = ('spearman r = %1.3f, p(spin) = %1.4f' % (corrval, pval))
xlab = 'Timescale - HCPD'
ylab = 'Timescale - HBN'
# fcn_timescale.scatterregplot(x, y, title, xlab, ylab, 50)
myplot = sns.scatterplot(x=x, y=y, palette=cmap_YlPu,
                         hue=sa_rank['SArank'].values,
                         legend=True)
sns.regplot(x=x, y=y, scatter=False, ax=myplot,
            line_kws=dict(color='k'))
sns.despine(ax=myplot, trim=False)
myplot.axes.set_title(title)
myplot.axes.set_xlabel(xlab)
myplot.axes.set_ylabel(ylab)
myplot.figure.set_figwidth(6)
myplot.figure.set_figheight(6)
plt.tight_layout()
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_timescale_corr_v2.svg',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_timescale_corr_v2.png',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.close()

################
# higher order vs lower order
################
hbn_ts = pd.read_csv(inpath + 'data/timescale/' +
                     'HBN_rest_noSI_age_TRorig_concat/' +
                     'HBN_rest_noSI_age_concat_timescale_acf_TRorig' +
                     '_Schaefer_400-7_forR.tsv', sep='\t')
hcpd_ts = pd.read_csv(inpath + 'data/timescale/HCPD_rest_age_TRorig_concat/' +
                      'HCPD_rest_age_concat_timescale_acf_TRorig' +
                      '_Schaefer_400-7_forR.tsv', sep='\t')
hcpya_ts = pd.read_csv(inpath + 'data/timescale/HCPYA_rest_TRorig_concat/' +
                       'HCPYA_rest_concat_timescale_acf_TRorig' +
                       '_Schaefer_400-7_forR.tsv', sep='\t')


x = pd.DataFrame(np.nanmean(np.array(hcpd_ts.iloc[:, :400]), axis=0))
y = pd.DataFrame(np.nanmean(np.array(hbn_ts.iloc[:, :400]), axis=0))
z = pd.DataFrame(np.nanmean(np.array(hcpya_ts.iloc[:, :400]), axis=0))

combined_df = pd.concat([x, y, z], axis=1)
labels = ['hcpd', 'hbn', 'hcpya']
combined_df.columns = labels

combined_df['SArank'] = sa_rank['SArank'].values
combined_df['category'] = combined_df['SArank'].apply(lambda x: 'sensorimotor'
                                                      if 1 <= x <= 200
                                                      else 'association')

dataset = 'hcpya'
# Separate regions based on the SA category column
sensorimotor_regions = combined_df[combined_df['category'] ==
                                   'sensorimotor'][dataset]
association_regions = combined_df[combined_df['category'] ==
                                  'association'][dataset]

# Perform an independent t-test (Welch's t-test by default)
t_stat, p_value = ttest_ind(association_regions,
                            sensorimotor_regions, equal_var=False)

nspin = 10000
spin_idx = fcn_timescale.get_spinidx(nspin=nspin, lhannot=lhlabels,
                                     rhannot=rhlabels, surfpath=surf_path)
permuted_sa_rank = sa_rank['SArank'].values[spin_idx]
combined_df_perm = combined_df.copy()
t_stat_null = []

for sIdx in np.arange(nspin):
    combined_df_perm['SArank_perm'] = permuted_sa_rank[:, sIdx]
    combined_df_perm['category_perm'] = combined_df_perm['SArank_perm'].apply(
        lambda x: 'sensorimotor' if 1 <= x <= 200 else 'association')
    sm_regions = combined_df_perm[combined_df_perm['category_perm'] ==
                                  'sensorimotor'][dataset]
    assoc_regions = combined_df_perm[combined_df_perm['category_perm'] ==
                                     'association'][dataset]
    t_stat_perm, _ = ttest_ind(assoc_regions, sm_regions, equal_var=False)
    t_stat_null.append(t_stat_perm)

permmean = np.mean(t_stat_null)
pvalspin = (len(np.where(abs(t_stat_null - permmean) >=
                         abs(t_stat - permmean))[0])+1)/(nspin+1)

to_plot = pd.concat([sensorimotor_regions, association_regions], axis=1)
labels = ['sensorimotor', 'association']
to_plot.columns = labels

plt.ion()
sns.set_theme(style='ticks', palette='pastel')
flierprops = dict(marker='.')
myplot = sns.boxplot(data=to_plot, palette=['m', 'g'],
                     width=0.4, fliersize=3, showcaps=False,
                     flierprops=flierprops, showfliers=False)
sns.despine(ax=myplot, offset=5, trim=True)
myplot.axes.set_title('%s: t-value = %1.3f, - p (spin) = %1.4f'
                      % (dataset.upper(), t_stat, pvalspin))
myplot.set(ylabel='timescale')
myplot.figure.set_figwidth(4)
myplot.figure.set_figheight(4)
plt.tight_layout()
plt.show()
plt.savefig(inpath + 'results/timescale/' +
            'timescale_SA_comparison_%s.svg' % dataset,
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.savefig(inpath + 'results/timescale/' +
            'timescale_SA_comparison_%s.png' % dataset,
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.close()

################
# correlate TS with SA rank with spins
################
hbn_ts = pd.read_csv(inpath + 'data/timescale/' +
                     'HBN_rest_noSI_age_TRorig_concat/' +
                     'HBN_rest_noSI_age_concat_timescale_acf_TRorig' +
                     '_Schaefer_400-7_forR.tsv', sep='\t')
hcpd_ts = pd.read_csv(inpath + 'data/timescale/HCPD_rest_age_TRorig_concat/' +
                      'HCPD_rest_age_concat_timescale_acf_TRorig' +
                      '_Schaefer_400-7_forR.tsv', sep='\t')
hcpya_ts = pd.read_csv(inpath + 'data/timescale/HCPYA_rest_TRorig_concat/' +
                       'HCPYA_rest_concat_timescale_acf_TRorig' +
                       '_Schaefer_400-7_forR.tsv', sep='\t')


x = pd.DataFrame(np.nanmean(np.array(hcpd_ts.iloc[:, :400]), axis=0))
y = pd.DataFrame(np.nanmean(np.array(hbn_ts.iloc[:, :400]), axis=0))
z = pd.DataFrame(np.nanmean(np.array(hcpya_ts.iloc[:, :400]), axis=0))

combined_df = pd.concat([x, y, z], axis=1)
labels = ['hcpd', 'hbn', 'hcpya']
combined_df.columns = labels

combined_df['SArank'] = sa_rank['SArank'].values

# spin p
nspins = 10000
corrval = []
pvalspin = []
for iVar in range(len(labels)):
    x = sa_rank['SArank'].values
    y = combined_df.iloc[:, iVar].values
    corr_temp = scipy.stats.spearmanr(x, y)[0]
    pval_temp = fcn_timescale.get_spinp(x, y, corrval=corr_temp, nspin=nspins,
                                        lhannot=lhlabels, rhannot=rhlabels,
                                        corrtype='spearman',
                                        surfpath=surf_path)
    corrval.append(corr_temp)
    pvalspin.append(pval_temp)

sa_spin = pd.DataFrame({'corrval': corrval,
                        'pvalspin': pvalspin,
                        'dataset': labels})
sa_spin.to_csv(inpath + 'results/timescale/csvFiles/sa_ts_spin_results.csv',
               index=False)

################
# correlate R2 from R
################
hbn = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                  'HBN_rest_noSI_age_timescale_age_r2_TRorig.csv')
hcpd = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                   'HCPD_rest_age_timescale_age_r2_TRorig.csv')

combined_df = pd.concat([hcpd['partialR2'],
                         hbn['partialR2']], axis=1)
labels = ['hcpd', 'hbn']
combined_df.columns = labels
combined_data = np.array(combined_df)

# copy region info
region_names = regionlabels.copy()

region_names_r2 = list(hcpd['gam.age.schaefer.region'].values)

sort_idx = []
for region in region_names:
    temp_idx = np.where(np.array(region_names_r2) == region)[0][0]
    sort_idx.append(temp_idx)

combined_data_reordered = combined_data.copy()
combined_data_reordered = combined_data_reordered[sort_idx, :]

combined_data_reordered = pd.DataFrame(combined_data_reordered,
                                       columns=labels)

combined_data_reordered['hcpd-rank'] = combined_data_reordered['hcpd'].rank()
combined_data_reordered['hbn-rank'] = combined_data_reordered['hbn'].rank()

# spin p
nspins = 10000

x = combined_data_reordered['hcpd'].values
y = combined_data_reordered['hbn'].values

corrval = scipy.stats.spearmanr(x, y)[0]
pval = fcn_timescale.get_spinp(x, y, corrval=corrval, nspin=nspins,
                               lhannot=lhlabels, rhannot=rhlabels,
                               corrtype='spearman',
                               surfpath=surf_path)
plt.ion()
title = ('spearman r = %1.3f, p(spin) = %1.4f' % (corrval, pval))
xlab = 'partial R2 - HCPD'
ylab = 'partial R2 - HBN'
# fcn_timescale.scatterregplot(x, y, title, xlab, ylab, 50)
myplot = sns.scatterplot(x=x, y=y, palette=cmap_YlPu,
                         hue=sa_rank['SArank'].values,
                         legend=True)
sns.regplot(x=x, y=y, scatter=False, ax=myplot,
            line_kws=dict(color='k'))
sns.despine(ax=myplot, trim=False)
myplot.axes.set_title(title)
myplot.axes.set_xlabel(xlab)
myplot.axes.set_ylabel(ylab)
myplot.figure.set_figwidth(6)
myplot.figure.set_figheight(6)
plt.tight_layout()
plt.tight_layout()
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_partialR2_corr_v2.svg',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_partialR2_corr_v2.png',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.close()

################
# correlate R2 from R with SA rank with spins
################
hbn = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                  'HBN_rest_noSI_age_timescale_age_r2_TRorig.csv')
hcpd = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                   'HCPD_rest_age_timescale_age_r2_TRorig.csv')
hcpya = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                    'HCPYA_rest_timescale_age_r2_TRorig.csv')

combined_df = pd.concat([hbn['partialR2'],
                         hcpd['partialR2'],
                         hcpya['partialR2']], axis=1)
labels = ['hbn', 'hcpd', 'hcpya']
combined_data = np.array(combined_df)

# copy region info
region_names = regionlabels.copy()

region_names_r2 = list(hbn['gam.age.schaefer.region'].values)

sort_idx = []
for region in region_names:
    temp_idx = np.where(np.array(region_names_r2) == region)[0][0]
    sort_idx.append(temp_idx)

combined_data_reordered = combined_data.copy()
combined_data_reordered = combined_data_reordered[sort_idx, :]

# spin p
nspins = 10000
corrval = []
pvalspin = []
for iVar in range(len(labels)):
    x = sa_rank['SArank'].values
    y = combined_data_reordered[:, iVar]
    corr_temp = scipy.stats.spearmanr(x, y)[0]
    pval_temp = fcn_timescale.get_spinp(x, y, corrval=corr_temp, nspin=nspins,
                                        lhannot=lhlabels, rhannot=rhlabels,
                                        corrtype='spearman',
                                        surfpath=surf_path)
    corrval.append(corr_temp)
    pvalspin.append(pval_temp)

sa_spin = pd.DataFrame({'corrval': corrval,
                        'pvalspin': pvalspin,
                        'dataset': labels})
sa_spin.to_csv(inpath + 'results/timescale/csvFiles/sa_r2_spin_results.csv',
               index=False)

###################
# brain volume
###################
hbn_ts = pd.read_csv(inpath + 'data/timescale/' +
                     'HBN_rest_noSI_age_TRorig_concat/' +
                     'HBN_rest_noSI_age_concat_timescale_acf_TRorig' +
                     '_Schaefer_400-7_forR.tsv', sep='\t')
hcpd_ts = pd.read_csv(inpath + 'data/timescale/HCPD_rest_age_TRorig_concat/' +
                      'HCPD_rest_age_concat_timescale_acf_TRorig' +
                      '_Schaefer_400-7_forR.tsv', sep='\t')

hcpd_vol = pd.read_csv(inpath + 'data/brain_volume/df_vol_hcpd.csv')

x = pd.DataFrame(np.nanmean(np.array(hbn_ts.iloc[:, :400]), axis=0))
y = pd.DataFrame(np.nanmean(np.array(hcpd_ts.iloc[:, :400]), axis=0))

subjList_ts = [int(i.split('-')[1]) for i in list(hcpd_ts['participant_id'])]
subjList_vol = hcpd_vol['subjID']

x = np.array(subjList_ts)
y = np.array(subjList_vol)
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

ts_shared = hcpd_ts.iloc[x_ind, :]
ts_shared.reset_index(inplace=True, drop=True)

vol_shared = hcpd_vol.iloc[y_ind, :]
vol_shared.reset_index(inplace=True, drop=True)

x = np.nanmean(np.array(ts_shared.iloc[:, :400]), axis=1)
y = np.array(vol_shared['vol'])

corrvol = scipy.stats.spearmanr(x, y)[0]
print('\nSpearman r between total brain volume and timescale in HCPD: %.4f'
      % corrvol)

# calculate parcel size
lh_data = nib.load(lhlabels)
rh_data = nib.load(rhlabels)

lh_data = lh_data.darrays[0].data
rh_data = rh_data.darrays[0].data

parcelsize_lh = [np.count_nonzero(lh_data == int(i))
                 for i in np.unique(lh_data)[1:]]
parcelsize_rh = [np.count_nonzero(rh_data == int(i))
                 for i in np.unique(rh_data)[1:]]
parcelsize = parcelsize_lh + parcelsize_rh
parcelsize = np.array(parcelsize)

timescale_data = np.array(hcpd_ts.iloc[:, :400])
hcpd_corr = []
for iRow in range(len(hcpd_ts)):
    subj_ts = timescale_data[iRow, :]
    hcpd_corr.append(scipy.stats.pearsonr(parcelsize, subj_ts)[0])

timescale_data = np.array(hbn_ts.iloc[:, :400])
hbn_corr = []
for iRow in range(len(hcpd_ts)):
    subj_ts = timescale_data[iRow, :]
    hbn_corr.append(scipy.stats.pearsonr(parcelsize, subj_ts)[0])

print('\nSpearman r mean parcel size and timescale in HCPD: %.4f'
      % np.mean(hcpd_corr))
print('\nSpearman r mean parcel size and timescale in HBN: %.4f'
      % np.mean(hbn_corr))

# age effect
hbn = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                  'HBN_rest_noSI_age_timescale_age_r2_TRorig.csv')
hcpd = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                   'HCPD_rest_age_timescale_age_r2_TRorig.csv')

# copy region info
region_names = regionlabels.copy()
region_names_r2 = list(hbn['gam.age.schaefer.region'].values)

sort_idx = []
for region in region_names:
    temp_idx = np.where(np.array(region_names_r2) == region)[0][0]
    sort_idx.append(temp_idx)

ageeffect_data = np.array(hcpd['partialR2'])
ageeffect_data_reordered = ageeffect_data[sort_idx]
hcpd_corr_ageeffect = scipy.stats.pearsonr(parcelsize,
                                           ageeffect_data_reordered)[0]

ageeffect_data = np.array(hbn['partialR2'])
ageeffect_data_reordered = ageeffect_data[sort_idx]
hbn_corr_ageeffect = scipy.stats.pearsonr(parcelsize,
                                          ageeffect_data_reordered)[0]

print('\nSpearman r mean parcel size and age effect (R2) in HCPD: %.4f'
      % np.mean(hcpd_corr_ageeffect))
print('\nSpearman r mean parcel size and age effect (R2) in HBN: %.4f'
      % np.mean(hbn_corr_ageeffect))

################
# correlate other age stats from R -- will prob go
################
hbn = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                  'HBN_rest_noSI_age_timescale_age_statistics_full_TRorig.csv')
hcpd = pd.read_csv(inpath + 'results/timescale/csvFiles/' +
                   'HCPD_rest_age_timescale_age_statistics_full_TRorig.csv')

hcpd['region'] = hcpd['region'].str.removeprefix('X')
hbn['region'] = hbn['region'].str.removeprefix('X')

# age stat of interest: GAM.age.rankR2sig, age.maturation, age.peakchange
combined_df = pd.concat([hcpd['GAM.age.rankR2sig'],
                         hbn['GAM.age.rankR2sig']], axis=1)
labels = ['hcpd', 'hbn']
combined_df.columns = labels
combined_data = np.array(combined_df)

# copy region info
region_names = regionlabels.copy()

region_names_r2 = list(hcpd['region'].values)

sort_idx = []
for region in region_names:
    temp_idx = np.where(np.array(region_names_r2) == region)[0][0]
    sort_idx.append(temp_idx)

combined_data_reordered = combined_data.copy()
combined_data_reordered = combined_data_reordered[sort_idx, :]

combined_data_reordered = pd.DataFrame(combined_data_reordered,
                                       columns=labels)

x = combined_data_reordered['hcpd'].values
y = combined_data_reordered['hbn'].values


def positional_overlap_score(x, y, mode='jaccard'):
    valid_x = ~np.isnan(x)
    valid_y = ~np.isnan(y)
    intersection = np.sum(valid_x & valid_y)
    union = np.sum(valid_x | valid_y)
    if mode == 'jaccard':
        return intersection / union
    elif mode == 'dice':
        return 2 * intersection / (np.sum(valid_x) + np.sum(valid_y))


# spatial overlap
overlap_score = positional_overlap_score(x, y, mode='dice')
print('overlap score:', overlap_score)

# value agreement
mask = ~np.isnan(x) & ~np.isnan(y)
corrval = scipy.stats.spearmanr(x[mask], y[mask])[0]
print('spearman r:', corrval)

# composite score for spatial overlap and value agreement
composite_score = corrval * overlap_score
print('composite score:', composite_score)

nspin = 10000
spin_idx = fcn_timescale.get_spinidx(nspin=nspin, lhannot=lhlabels,
                                     rhannot=rhlabels, surfpath=surf_path)
permuted_x = x[spin_idx]
overlap_score_null = []
corrval_null = []
composite_score_null = []

for sIdx in np.arange(nspin):
    # spatial overlap
    overlap_score_perm = positional_overlap_score(permuted_x[:, sIdx], y,
                                                  mode='dice')
    # value agreement
    mask = ~np.isnan(permuted_x[:, sIdx]) & ~np.isnan(y)
    corrval_perm = scipy.stats.spearmanr(permuted_x[:, sIdx][mask], y[mask])[0]

    composite_score_perm = corrval_perm * overlap_score_perm

    overlap_score_null.append(overlap_score_perm)
    corrval_null.append(corrval_perm)
    composite_score_null.append(composite_score_perm)

permmean = np.mean(overlap_score_null)
pvalspin_overlap = (len(np.where(abs(overlap_score_null - permmean) >=
                                 abs(overlap_score - permmean))[0])+1)/(nspin
                                                                        + 1)

permmean = np.mean(corrval_null)
pvalspin_corrval = (len(np.where(abs(corrval_null - permmean) >=
                                 abs(corrval - permmean))[0])+1)/(nspin+1)

permmean = np.mean(composite_score_null)
pvalspin_composite = (len(np.where(abs(composite_score_null - permmean) >=
                                   abs(composite_score - permmean))[0])+1)/(
                                       nspin+1)

print('overlap score:', overlap_score, ', pval(spin):', pvalspin_overlap)
print('spearman r:', corrval, ', pval(spin):', pvalspin_corrval)
print('composite score:', composite_score, ', pval(spin):', pvalspin_composite)

plt.ion()
ax = sns.displot(overlap_score_null, kind='kde')
# ax.set_ylim([0, 6])
plt.vlines(overlap_score, ymin=0, ymax=11, colors='r')
plt.xlabel('spatial overlap (Dice index)')
plt.ylabel('kernel density estimate')
plt.title('FDR corrected partial R2: Dice index = %1.3f, p(spin) = %1.4f'
          % (overlap_score, pvalspin_overlap))
plt.tight_layout()
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_fdrR2_dice.svg',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_fdrR2_dice.png',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.close()

# or just correlate and plot
nspins = 10000
x = combined_data_reordered['hcpd'].values
y = combined_data_reordered['hbn'].values

corrval = scipy.stats.spearmanr(x, y)[0]
pval = fcn_timescale.get_spinp(x, y, corrval=corrval, nspin=nspins,
                               lhannot=lhlabels, rhannot=rhlabels,
                               corrtype='spearman',
                               surfpath=surf_path)

plt.ion()
title = ('spearman r = %1.3f, p(spin) = %1.4f' % (corrval, pval))
title = ('spearman r = %1.3f' % corrval)
xlab = 'Age of maturation - HCPD'
ylab = 'Age of maturation - HBN'
fcn_timescale.scatterregplot(x, y, title, xlab, ylab, 50)
plt.tight_layout()
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_agematuration_corr.svg',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.savefig(inpath + 'results/timescale/' +
            'hbn_hcpd_agematuration_corr.png',
            bbox_inches='tight', dpi=300,
            transparent=True)
plt.close()
