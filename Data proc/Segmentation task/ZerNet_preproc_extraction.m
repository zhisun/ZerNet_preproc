function ZerNet_preproc_extraction(srcroot, dstroot_1, dstroot_2, scale_infor)
filelist = dir(fullfile(srcroot, '*.mat'));
if ~exist(dstroot_1, 'dir')
    mkdir(dstroot_1);
end
if ~exist(dstroot_2, 'dir')
    mkdir(dstroot_2);
end

parfor model_id =1:length(filelist)
    
    fprintf('[I] ZerNet_preproc extraction for model %s\n\n', filelist(model_id).name);
    pro_start_time = tic;
    
    mesh_label=load(fullfile(filelist(model_id).folder, filelist(model_id).name));
    shape = struct('X',mesh_label.X(:,1),'Y',mesh_label.X(:,2),'Z',mesh_label.X(:,3),'TRIV',mesh_label.F,'label_V',mesh_label.label_X,'label_F',mesh_label.label_F);
    mesh_label = [];

    [ZerNet_preproc,sampled_surface] = extract_ZerPatch(scale_infor,shape);
    
    savepath_1 = fullfile(dstroot_1,filelist(model_id).name);
    parsave_1(savepath_1,sampled_surface);
    
    savepath_2 = fullfile(dstroot_2,filelist(model_id).name);
    parsave_2(savepath_2,ZerNet_preproc);
    
    elapsed_time = toc(pro_start_time);
    fprintf('model %s ZerNet_preproc extraction time: %2.4fs\n',filelist(model_id).name, elapsed_time);
end
end

function parsave_1(path,sampled_surface)
save(path,'sampled_surface')
end

function parsave_2(path,ZerNet_preproc)
save(path,'ZerNet_preproc')
end