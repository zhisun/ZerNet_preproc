function [Input_ZerCoeff, Recons_acc] = compute_Input_Zercoeff(shape,resample_graph,disk_R,n_rays,zernike_order,num_bases,split_num,mean_disk_num)
if split_num < num_bases
    error('Error. \split_num should be equal or greater than num of Zernike-bases.')
end
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
N_X = resample_graph.N_X;
D_X2N = resample_graph.D_X2N;
%
% [N_X,I,B,r] = UniformSampling_on_mesh(X,F,resample_graph.sampled_num,'Color','blue');
% resample_graph.N_X = N_X;
% resample_graph.I = I;
% resample_graph.B = B;
% resample_graph.r = r;
% %compute geodiesic distance between original mesh vertices and sampled surface points
% D_X2N = X2N_GeodesicDistance(F,D,I,B);
%
N = []; M = [];
for n = 0:zernike_order
    N = [N n*ones(1,n+1)];
    M = [M -n:2:n];
end

% X_Disk_cors = cell(length(X),1);
% X_Disk_dis = cell(length(X),1);
% X_Disk_theta = cell(length(X),1);
% X_Disk_Zerbases = cell(length(X),1);
% X_Disk_FeaVec = cell(length(X),1);

X_Disk_zernike_coeff = zeros(size(X,1),num_bases,size(X,2));
X_Disk_recons_acc = zeros(size(X,1),1);
%
vfring = compute_vertex_face_ring(F');
shape.idxs = vfring;
shape.f_dns = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));

for Xid = 1:length(X)
    center_cors = X(Xid,:)';
    D_Xid2N = D_X2N(Xid,:);
    Disk_Nids = find(D_Xid2N < disk_R & D_Xid2N > 0);
    Disk_Ndis =  D_Xid2N(Disk_Nids);
    %
    [direction, theta] = get_sample_directions(X', F', vfring, Xid, n_rays);
    shape = fast_marching(Xid,shape,'vertex',0,1,shape.f_dns);
    
    Disk_Xids = find(shape.D < disk_R & shape.D > 0);
    Disk_Xdis =  shape.D(Disk_Xids);
    %    
    cor_intersections = [];
    theta_intersection = [];
    for i=1:n_rays
        [~, ray_intersections] = geodesic_triangles(Xid, direction(:,i), shape, 0.1);%disk_R);
        ray_intersections = ray_intersections(:,2:end);
        ray_thetas = repmat(theta(i), 1,size(ray_intersections,2));
        cor_intersections = [cor_intersections,ray_intersections];
        theta_intersection = [theta_intersection,ray_thetas];
    end
    vec_intersections = cor_intersections - repmat(center_cors,1,length(cor_intersections));
    if length(Disk_Nids)<split_num
        filled_num = min((mean_disk_num - length(Disk_Nids)),length(Disk_Xdis));
        pick_ids = randsample(length(Disk_Xdis),filled_num);
        Disk_Xids =  Disk_Xids(pick_ids);
        Disk_Xdis = Disk_Xdis(pick_ids);
        clear pick_ids
        Disk_cors = [N_X(Disk_Nids,:)',X(Disk_Xids,:)'];
        Disk_dis = [Disk_Ndis,Disk_Xdis'];
    else
        Disk_cors = N_X(Disk_Nids,:)';
        Disk_dis = Disk_Ndis;
    end
    [Disk_dis,tem_idx] = sort(Disk_dis);
    Disk_cors = Disk_cors(:,tem_idx);
    clear tem_idx
    
    vec_Disk_cors = Disk_cors - repmat(center_cors,1,length(Disk_cors)); 
    [~,Idxs]=pdist2(vec_intersections',vec_Disk_cors','cosine','Smallest',1);
    Disk_theta = theta_intersection(Idxs);
    
%     X_Disk_cors{Xid} = Disk_cors;
%     X_Disk_dis{Xid} = Disk_dis;
%     X_Disk_theta{Xid} = Disk_theta;

%     figure; plot_mesh(X, F);hold on;
%     scatter3(X(Xid,1),X(Xid,2),X(Xid,3), 20, 'filled', 'r')
%     vis_ids_1 = find(Disk_theta>=0 & Disk_theta<pi);
%     vis_ids_2 = find(Disk_theta>=pi & Disk_theta<2*pi);
%     scatter3(Disk_cors(1,vis_ids_1),Disk_cors(2,vis_ids_1),Disk_cors(3,vis_ids_1), 2, 'filled', 'g');
%     scatter3(Disk_cors(1,vis_ids_2),Disk_cors(2,vis_ids_2),Disk_cors(3,vis_ids_2), 2, 'filled', 'y');
    
    %calculate Zernike bases in each disk
    r_k = Disk_dis/disk_R;
    theta_k = Disk_theta;
    Disk_bases = zernfun(N,M,r_k,theta_k,'norm');
    %calculate input feature vec in each disk
    Disk_FeaVec = vec_Disk_cors;
    Disk_zernike_coeff = Disk_bases\Disk_FeaVec';
    
    recon_Disk_FeaVec = Disk_bases*Disk_zernike_coeff;
    recon_error =  abs((recon_Disk_FeaVec - Disk_FeaVec')./Disk_FeaVec');
    recons_acc = length(find(recon_error < 0.2))/length(recon_error(:));
    
%     X_Disk_Zerbases{Xid} = Disk_bases;
%     X_Disk_FeaVec{Xid} = Disk_FeaVec';
    
    X_Disk_zernike_coeff(Xid,:,:) = Disk_zernike_coeff;
    X_Disk_recons_acc(Xid) = recons_acc;
end
Input_ZerCoeff = X_Disk_zernike_coeff;
Recons_acc = X_Disk_recons_acc;
end