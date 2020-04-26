%% input data preprocess for Faust dataset for point-wise correspondence task
% Zernike patches will be computed.
mesh_root = '../../../ZerNet/Data/Faust Corrspondence Dataset/meshes/area=15000';
%
mapping_infor_root = '../../../ZerNet/Data/Faust Corrspondence Dataset/ZerNet Input/X2N_mapping_infor';
Zer_patches_root = '../../../ZerNet/Data/Faust Corrspondence Dataset/ZerNet Input/Input ZerPatches';
scale_infor.ref_model_id= 1;

% single_scale
scale_infor.scale_1.ref_sample_num = 12000;
scale_infor.scale_1.disk_R = 5.5;
scale_infor.scale_1.n_rays = 36;
scale_infor.scale_1.cut_num = 50;
scale_infor.scale_1.zernike_order = 5;
scale_infor.scale_1.num_bases = 21;

%-------------------------------------------------------------------------------------------------------
% % multi_sacle: Used for multiscale ZerNet, at three different size of
% receptive field, Zernike patches at each scale will be computed.

% mapping_infor_root = '../../../Data/Faust Dataset/ZerNet Input/X2N_mapping_infor/R=4.5_5.5_6.75';
% Zer_patches_root = '../../../Data/Faust Dataset/ZerNet Input/Input ZerPatches/R=4.5_5.5_6.75';
% scale_infor.ref_model_id= 1;

% % scale_1
% scale_infor.scale_1.ref_sample_num =18000;
% scale_infor.scale_1.disk_R = 4.5;
% scale_infor.scale_1.n_rays = 36;
% scale_infor.scale_1.cut_num = 50;
% scale_infor.scale_1.zernike_order = 5;
% scale_infor.scale_1.num_bases = 21;
% 
% % scale_2
% scale_infor.scale_2.ref_sample_num = 12000;
% scale_infor.scale_2.disk_R = 5.5;
% scale_infor.scale_2.n_rays = 36;
% scale_infor.scale_2.cut_num = 50;
% scale_infor.scale_2.zernike_order = 5;
% scale_infor.scale_2.num_bases = 21;
% 
% % scale_3
% scale_infor.scale_3.ref_sample_num = 8000;
% scale_infor.scale_3.disk_R = 6.75;
% scale_infor.scale_3.n_rays = 36;
% scale_infor.scale_3.cut_num = 50;
% scale_infor.scale_3.zernike_order = 5;
% scale_infor.scale_3.num_bases = 21;
%------------------------------------------------------------------------------------------------------------------
ZerNet_preproc_extraction(mesh_root, mapping_infor_root, Zer_patches_root,scale_infor);
%% Extract labels as ZerNet Outputs
Output_root = '../../../ZerNet/Data/Faust Corrspondence Dataset/ZerNet Output';
num_Vs = 6890;
extract_label(mesh_root, Output_root,num_Vs);

%% Extract mesh vertices XYZ coordinates as ZerNet initial input features
feature_input_root = '../../../ZerNet/Data/Faust Corrspondence Dataset/ZerNet Input/Input features';
extract_feature(mesh_root, feature_input_root);
