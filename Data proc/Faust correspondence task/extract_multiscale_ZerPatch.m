function ZerNet_preproc = extract_multiscale_ZerPatch(scale_infor,shape)
fprintf('Compute geodesic distance of mesh vertices, elapsed_time: ');
start_time  = tic;
shape.D_global = compute_geods(shape);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

% scale_1
scale_params = scale_infor.scale_1;
ref_resample_graph = scale_params.ref_resample_graph;
disk_R = scale_params.disk_R;
n_rays = scale_params.n_rays; 
cut_num = scale_params.cut_num;
zernike_order = scale_params.zernike_order;
num_bases = scale_params.num_bases;
clear scale_params

%
fprintf('Compute scale_1 resample graph via precomputed interpolation weights on original mesh vertices, elapsed_time: ');
start_time  = tic;
resample_graph = construct_resample_graph_via_ref(shape,ref_resample_graph);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

fprintf('Compute Zernike patches at scale_1, elapsed_time: ');
start_time  = tic;
ZerPatch = compute_ZerPatch(shape,resample_graph,disk_R,n_rays,cut_num,zernike_order,num_bases);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

ZerNet_preproc.scale_1.numofdPts = resample_graph.NumPids;
ZerNet_preproc.scale_1.patch_radius = disk_R;
ZerNet_preproc.scale_1.patch_numofdPts = cut_num;
ZerNet_preproc.scale_1.zernike_order = zernike_order;
ZerNet_preproc.scale_1.ZerPatch = ZerPatch;

clear scale_params ref_resample_graph resample_graph disk_R n_rays cut_num zernike_order num_bases ZerPatch

% %used for multiscle
% % scale_2
% scale_params = scale_infor.scale_2;
% ref_resample_graph = scale_params.ref_resample_graph;
% disk_R = scale_params.disk_R;
% n_rays = scale_params.n_rays; 
% cut_num = scale_params.cut_num;
% zernike_order = scale_params.zernike_order;
% num_bases = scale_params.num_bases;
% clear scale_params
% 
% %
% fprintf('Compute scale_2 resample graph via precomputed interpolation weights on original mesh vertices, elapsed_time: ');
% start_time  = tic;
% resample_graph = construct_resample_graph_via_ref(shape,ref_resample_graph);
% elapsed_time = toc(start_time);
% fprintf('%2.4fs\n',elapsed_time);
% 
% fprintf('Compute Zernike patches at scale_2, elapsed_time: ');
% start_time  = tic;
% ZerPatch = compute_ZerPatch(shape,resample_graph,disk_R,n_rays,cut_num,zernike_order,num_bases);
% elapsed_time = toc(start_time);
% fprintf('%2.4fs\n',elapsed_time);
% 
% ZerNet_preproc.scale_2.numofdPts = resample_graph.NumPids;
% ZerNet_preproc.scale_2.patch_radius = disk_R;
% ZerNet_preproc.scale_2.patch_numofdPts = cut_num;
% ZerNet_preproc.scale_2.zernike_order = zernike_order;
% ZerNet_preproc.scale_2.ZerPatch = ZerPatch;
% 
% clear scale_params ref_resample_graph resample_graph disk_R n_rays cut_num zernike_order num_bases ZerPatch
% 
% % scale_3
% scale_params = scale_infor.scale_3;
% ref_resample_graph = scale_params.ref_resample_graph;
% disk_R = scale_params.disk_R;
% n_rays = scale_params.n_rays; 
% cut_num = scale_params.cut_num;
% zernike_order = scale_params.zernike_order;
% num_bases = scale_params.num_bases;
% clear scale_params
% 
% %
% fprintf('Compute scale_3 resample graph via precomputed interpolation weights on original mesh vertices, elapsed_time: ');
% start_time  = tic;
% resample_graph = construct_resample_graph_via_ref(shape,ref_resample_graph);
% elapsed_time = toc(start_time);
% fprintf('%2.4fs\n',elapsed_time);
% 
% fprintf('Compute Zernike patches at scale_3, elapsed_time: ');
% start_time  = tic;
% ZerPatch = compute_ZerPatch(shape,resample_graph,disk_R,n_rays,cut_num,zernike_order,num_bases);
% elapsed_time = toc(start_time);
% fprintf('%2.4fs\n',elapsed_time);
% 
% ZerNet_preproc.scale_3.numofdPts = resample_graph.NumPids;
% ZerNet_preproc.scale_3.patch_radius = disk_R;
% ZerNet_preproc.scale_3.patch_numofdPts = cut_num;
% ZerNet_preproc.scale_3.zernike_order = zernike_order;
% ZerNet_preproc.scale_3.ZerPatch = ZerPatch;
% 
% clear scale_params ref_resample_graph resample_graph disk_R n_rays cut_num zernike_order num_bases ZerPatch

end