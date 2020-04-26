%% Input data preprocess for Aneurysm dataset for wall stress estimation task
% The wall stress at each aneurysm mesh vetices is named as "J2m" in "sampled_suface", and will be served as ZerNet output.

rawdata_root = '../../../ZerNet/Data/Aneurysm Dataset/mesh_stress_data/area=100';
sampled_surfaces_root = '../../../ZerNet/Data/Aneurysm Dataset/UniformSampling_surfaces';
Zer_patches_root = '../../../ZerNet/Data/Aneurysm Dataset/ZerNet Input/Input ZerPatches';

scale_infor.sample_num = 8000;
% scale_1
scale_infor.scale_1.disk_R = 0.6;
scale_infor.scale_1.n_rays = 32;
scale_infor.scale_1.cut_num = 65;
scale_infor.scale_1.zernike_order = 5;
scale_infor.scale_1.num_bases = 21;
scale_infor.scale_1.min_dis = 0.15;

ZerNet_preproc_extraction(rawdata_root, sampled_surfaces_root, Zer_patches_root,scale_infor);