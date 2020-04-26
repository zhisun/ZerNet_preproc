function extract_feature(srcroot, dstroot)

if ~exist(dstroot, 'dir')
    mkdir(dstroot);
end

filelist = dir(fullfile(srcroot, '*.mat'));
for model_id =1:length(filelist)
    mesh=load(fullfile(filelist(model_id).folder, filelist(model_id).name));
    shape = mesh.shape;
    mesh = [];
    feature_in = [shape.X,shape.Y,shape.Z];
    savepath = fullfile(dstroot,filelist(model_id).name);
    save(savepath,'feature_in');
end
end