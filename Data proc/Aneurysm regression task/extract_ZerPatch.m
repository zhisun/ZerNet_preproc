function [ZerNet_preproc,sampled_surface] = extract_ZerPatch(scale_infor,shape)
fprintf('Compute geodesic distance of mesh vertices, elapsed_time: ');
start_time  = tic;
shape.D_global = compute_geods(shape);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

sample_num = scale_infor.sample_num;
%
fprintf('Compute resample graph via uniformsampling on mesh surface, elapsed_time: ');
start_time  = tic;
resample_graph = compute_resample_graph(shape,sample_num);
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);
%
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
J2m = shape.J2m;
sampled_surface = struct('X',X,'F',F,'J2m',J2m,'I',resample_graph.I,'B',resample_graph.B,'N_X',resample_graph.N_X,'N_J2m',resample_graph.N_J2m);
clear X F J2m
%

% scale_1 
scale_params = scale_infor.scale_1;
disk_R = scale_params.disk_R;
n_rays = scale_params.n_rays; 
cut_num = scale_params.cut_num;
zernike_order = scale_params.zernike_order;
num_bases = scale_params.num_bases;
min_dis = scale_params.min_dis;
clear scale_params
%
% N_Disk_Nids_num = compute_disk_num(resample_graph,disk_R);
%
fprintf('Compute Zernike patches at scale_1, elapsed_time: ');
start_time  = tic;
ZerPatch = compute_ZerPatch(shape,resample_graph,disk_R,n_rays,cut_num,zernike_order,num_bases,min_dis);
%---------------------------------------------------------------------------------------------------------------------------
% This can be used while the ZeNet will directly take the resampled mesh
% surface as the target. Note that afther getting the ZerNet prediction on
% sampled surface point, it can be mapped back to orignial mesh vertices
% based on resampling information in "resample_graph" as a post-process.

% ZerPatch = compute_ZerPatch_resampled(shape, resample_graph, disk_R, n_rays, cut_num, zernike_order, num_bases, min_dis);
%----------------------------------------------------------------------------------------------------------------------------
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);

ZerNet_preproc.scale_1.patch_radius = disk_R;
ZerNet_preproc.scale_1.patch_numofdPts = cut_num;
ZerNet_preproc.scale_1.zernike_order = zernike_order;
ZerNet_preproc.scale_1.ZerPatch = ZerPatch;
clear disk_R n_rays cut_num zernike_order num_bases min_dis ZerPatch

end