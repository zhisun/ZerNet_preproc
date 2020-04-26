%% Input data preprocess for segmentaion task (Faust dataset is used as an example here)
%
rawdata_root = '../../../ZerNet/Data/Faust Segmentation Dataset/mesh_label_data/area=15000';
sampled_surfaces_root = '../../../ZerNet/Data/Faust Segmentation Dataset/UniformSampling_surfaces';
Zer_patches_root = '../../../ZerNet/Data/Faust Segmentation Dataset/ZerNet Input/Input ZerPatches';

scale_infor.sample_num = 12000;
% scale_1
scale_infor.scale_1.disk_R = 5.5;
scale_infor.scale_1.n_rays = 32;
scale_infor.scale_1.cut_num = 50;
scale_infor.scale_1.zernike_order = 5;
scale_infor.scale_1.num_bases = 21;
scale_infor.scale_1.min_dis = 0.5;

ZerNet_preproc_extraction(rawdata_root,sampled_surfaces_root,Zer_patches_root,scale_infor);