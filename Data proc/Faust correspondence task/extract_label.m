function extract_label(srcroot, dstroot,num_Vs)

if ~exist(dstroot, 'dir')
    mkdir(dstroot);
end

filelist = dir(fullfile(srcroot, '*.mat'));
for model_id =1:length(filelist)
    labels = (1:num_Vs)';
    savepath = fullfile(dstroot,filelist(model_id).name);
    save(savepath,'labels');
end
end