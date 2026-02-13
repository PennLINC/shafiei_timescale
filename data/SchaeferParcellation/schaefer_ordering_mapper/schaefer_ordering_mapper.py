import os
import wget
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
import scipy.io as sio

# %%
# params. note, nodes will be mapped from input_order to output_order
n_parcels = 200
input_order = 17
output_order = 7

# %% download data from Yeo lab github
# output dir
outdir = '/Users/gshafiei/Desktop/bpd/schaefer_ordering_mapper'

# github link
remote_path = 'https://github.com/ThomasYeoLab/CBIG/raw/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI'

# download schaefer files
for order in [input_order, output_order]:
    # nifti file in MNI space
    file = 'Schaefer2018_{0}Parcels_{1}Networks_order_FSLMNI152_2mm.nii.gz'.format(n_parcels, order)
    if os.path.exists(os.path.join(outdir, file)) == False:
        wget.download(os.path.join(remote_path, file), outdir)

    # centroids and labels
    file = 'Schaefer2018_{0}Parcels_{1}Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'.format(n_parcels, order)
    if os.path.exists(os.path.join(outdir, file)) == False:
        wget.download(os.path.join(remote_path, 'Centroid_coordinates', file), outdir)

# %% load centroids
file = 'Schaefer2018_{0}Parcels_{1}Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'.format(n_parcels, input_order)
centroids_inorder = pd.read_csv(os.path.join(outdir, file), index_col=0)

file = 'Schaefer2018_{0}Parcels_{1}Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'.format(n_parcels, output_order)
centroids_outorder = pd.read_csv(os.path.join(outdir, file), index_col=0)

# %% remap data
# mapped = pd.DataFrame(index=centroids_inorder.index, columns=centroids_inorder.columns)
mapped = pd.DataFrame(index=centroids_inorder.index)
for i, data in centroids_outorder.iterrows():
    # get coords to be matched from output order
    coords = [data['R'], data['A'], data['S']]

    # find index of matching node from input order
    idx = np.where(np.all(coords == centroids_inorder[['R', 'A', 'S']].values, axis=1))[0][0]

    # store
    mapped.loc[i, 'input_roi'] = centroids_inorder.loc[idx + 1, 'ROI Name']
    mapped.loc[i, 'output_roi'] = centroids_outorder.loc[i, 'ROI Name']
    mapped.loc[i, 'mapped_indices'] = idx

mapped['mapped_indices'] = mapped['mapped_indices'].astype(int)

# save out
mapped.to_csv(os.path.join(outdir, 'Schaefer_{0}-{1}_mappedto_{0}-{2}.csv'.format(n_parcels, input_order, output_order)))

# %%
# check accuracy using connectome's from Yeo's team
hcpdir = os.path.join(outdir, 'HCP')

fc_7network = sio.loadmat(os.path.join(hcpdir, '100206_Schaefer_7network_FC.mat'))
fc_7network = fc_7network['corr_mat'].copy()
fc_7network = fc_7network[:n_parcels, :][:, :n_parcels]

fc_17network = sio.loadmat(os.path.join(hcpdir, '100206_Schaefer_17network_FC.mat'))
fc_17network = fc_17network['corr_mat'].copy()
fc_17network = fc_17network[:n_parcels, :][:, :n_parcels]

# matrix of data stored in 17-network ordering
fc_17to7 = fc_17network.copy()
# matrix permuted to 7-network ordering
fc_17to7 = fc_17to7[mapped['mapped_indices'].values, :][:, mapped['mapped_indices'].values]

print(np.all(fc_7network == fc_17to7))
