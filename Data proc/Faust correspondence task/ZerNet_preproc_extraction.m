function ZerNet_preproc_extraction(srcroot, dstroot_1, dstroot_2,scale_infor)
filelist = dir(fullfile(srcroot, '*.mat'));
if ~exist(dstroot_1, 'dir')
    mkdir(dstroot_1);
end
ref_model_id = scale_infor.ref_model_id;
refmesh=load(fullfile(filelist(ref_model_id).folder, filelist(ref_model_id).name));
refshape = refmesh.shape;
X2N_mapping_infor.F = refshape.TRIV;
clear refmesh

fprintf('Compute geodesic distance of mesh vertices for ref_model %s\n\n', filelist(ref_model_id).name);
start_time  = tic;
refshape.D_global = compute_geods(refshape);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

%uniformsampling on reference model "tr_reg_000", and then applying the
%the same interpolation weights to all the other models .
fprintf('Compute ref_resample_graph via uniformsampling on ref_model %s\n\n', filelist(ref_model_id).name);
start_time  = tic;

%scale_1
ref_sample_num = scale_infor.scale_1.ref_sample_num;
scale_infor.scale_1.ref_resample_graph = compute_resample_graph(refshape,ref_sample_num);
clear ref_sample_num
X2N_mapping_infor.scale_1.I = scale_infor.scale_1.ref_resample_graph.I;
X2N_mapping_infor.scale_1.B = scale_infor.scale_1.ref_resample_graph.B;

% %%used for multiscle
% %scale_2
% ref_sample_num = scale_infor.scale_2.ref_sample_num;
% scale_infor.scale_2.ref_resample_graph = compute_resample_graph(refshape,ref_sample_num);
% scale_infor.scale_2.ref_resample_graph.D_X2N = [];
% clear ref_sample_num
% X2N_mapping_infor.scale_2.I = scale_infor.scale_2.ref_resample_graph.I;
% X2N_mapping_infor.scale_2.B = scale_infor.scale_2.ref_resample_graph.B;
% %scale_3
% ref_sample_num = scale_infor.scale_3.ref_sample_num;
% scale_infor.scale_3.ref_resample_graph = compute_resample_graph(refshape,ref_sample_num);
% scale_infor.scale_3.ref_resample_graph.D_X2N = [];
% clear ref_sample_num
% X2N_mapping_infor.scale_3.I = scale_infor.scale_3.ref_resample_graph.I;
% X2N_mapping_infor.scale_3.B = scale_infor.scale_3.ref_resample_graph.B;

elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);
clear refshape

savepath = fullfile(dstroot_1,'X2N_mapping_infor.mat');
save(savepath,'X2N_mapping_infor');
clear savepath

if ~exist(dstroot_2, 'dir')
    mkdir(dstroot_2);
end

parfor model_id =1:length(filelist)
    
    fprintf('[I] ZerNet_preproc extraction for model %s\n\n', filelist(model_id).name);
    pro_start_time = tic;
    mesh=load(fullfile(filelist(model_id).folder, filelist(model_id).name));
    shape = mesh.shape;
    mesh = [];
    ZerNet_preproc = extract_multiscale_ZerPatch(scale_infor,shape)
    savepath = fullfile(dstroot_2,filelist(model_id).name);
    parsave(savepath,ZerNet_preproc);
    elapsed_time = toc(pro_start_time);
    fprintf('model %s ZerNet_preproc extraction time: %2.4fs\n',filelist(model_id).name, elapsed_time);
end
end

function parsave(path,ZerNet_preproc)
save(path,'ZerNet_preproc')
end