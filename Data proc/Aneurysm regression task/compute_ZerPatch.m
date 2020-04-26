function ZerPatch = compute_ZerPatch(shape, resample_graph, disk_R, n_rays, cut_num, zernike_order, num_bases, min_r)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
D_X2N = X2N_GeodesicDistance(F,shape.D_global,resample_graph.I,resample_graph.B);
N_X = resample_graph.N_X;
I = resample_graph.I;
% cloest_N2X = knnsearch(X,N);
%
vfring = compute_vertex_face_ring(F');
shape.idxs = vfring;
shape.f_dns = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));

%compute sampled surface points normal
[V_normal,F_normal] = compute_normal(X,F);
N_normal = F_normal(:,I);
%
Nn = []; Mm = [];
for n = 0:zernike_order
    Nn = [Nn n*ones(1,n+1)];
    Mm = [Mm -n:2:n];
end
%
X_Disk_Nids = zeros(length(X),cut_num);
X_Disk_Zerbases = zeros(length(X),cut_num, num_bases);
X_Disk_FeaVec = zeros(length(X),cut_num, 3);
%

for Xid = 1:length(X)
    
    center_cors = X(Xid,:)';
    D_Xid2N = D_X2N(Xid,:);
    
    center_normal = V_normal(:,Xid);
    [~,nearest_Nid] = min(D_Xid2N);
    
    ref_dir = N_X(nearest_Nid,:)'-center_cors;
    projected_ref_dir = ref_dir - center_normal*dot(ref_dir,center_normal);
    projected_ref_dir = projected_ref_dir/norm(projected_ref_dir);
    
    Xid_faces=vfring{Xid};
    same_face_Nids = [];
    for i=1:length(Xid_faces)
        same_face_Nids = [same_face_Nids; find(I==Xid_faces(i))];
    end
    same_face_Nids = same_face_Nids';
    Dis = vecnorm(N_X(same_face_Nids,:)'-repmat(center_cors,1,length(same_face_Nids)));
    D_Xid2N(same_face_Nids) = Dis;
    
%     lower_bound = max(Dis);
%     lower_bound = min(lower_bound, min_r);
    lower_bound = min_r;
    
    Neibor_Nids = find(D_Xid2N < disk_R & D_Xid2N > lower_bound);
    Neibor_Ndis = D_Xid2N(Neibor_Nids);
    Neibor_dirs = N_X(Neibor_Nids,:)'-repmat(center_cors,1,length(Neibor_Nids));
    
    Neibor_thetas = zeros(1,length(Neibor_Nids));
    for i = 1:length(Neibor_Nids)
        projected_Neibor_dir = Neibor_dirs(:,i) - center_normal*dot(Neibor_dirs(:,i),center_normal);
        r = vrrotvec(projected_ref_dir,projected_Neibor_dir);
        rot_axis = r(1:3)';
        if(round(dot(rot_axis,center_normal))>0)
            Neibor_thetas(i)=r(4);
        else
            Neibor_thetas(i)=2*pi-r(4);
        end
    end
    
    Disk_Nids = Neibor_Nids;
    Disk_Ndis = Neibor_Ndis;
    Disk_theta = Neibor_thetas;
    
    if (length(Disk_Nids) < cut_num)
        if (2*length(Disk_Nids)> cut_num)
            duplicate_num = cut_num-length(Disk_Nids);
            ids_add = randsample(length(Disk_Nids),duplicate_num);
            Disk_Nids_add = Disk_Nids(ids_add);
            Disk_Ndis_add = Disk_Ndis(ids_add);
            rand_scale = (1.0-0.9).*rand(size(Disk_Ndis_add,1),size(Disk_Ndis_add,2))+0.9;
            Disk_Ndis_add = Disk_Ndis_add.*rand_scale;
            tem = randi(2,length(ids_add),1);
            rotated_interval = 2*pi/n_rays;
            THETA_add = Disk_theta(ids_add) + ((repmat(-1,length(ids_add),1).^tem)*rotated_interval)';
            Disk_Nids = [Disk_Nids,Disk_Nids_add];
            Disk_Ndis =[Disk_Ndis,Disk_Ndis_add];
            Disk_theta = [Disk_theta,THETA_add];
        else
            duplicate_ratio =  ceil((cut_num)/length(Disk_Nids));
            Disk_Nids = repmat(Disk_Nids, 1, duplicate_ratio);
            Disk_Ndis = repmat(Disk_Ndis, 1, duplicate_ratio);
            rand_scale = (1.0-0.9).*rand(size(Disk_Ndis,1),size(Disk_Ndis,2))+0.9;
            Disk_Ndis = Disk_Ndis.*rand_scale;
            
            Disk_theta = repmat(Disk_theta, 1, duplicate_ratio);
            rotated_list = [-2*pi/n_rays,0,2*pi/n_rays];
            Disk_theta = Disk_theta + randsample(rotated_list,length(Disk_theta),true);
        end        
        idx = find(Disk_theta < 0);
        Disk_theta(idx) = 2*pi + Disk_theta(idx);
        clear idx
        idx = find(Disk_theta > 2*pi);
        Disk_theta(idx) = Disk_theta(idx) - 2*pi;
        clear idx
    end
    
    [Disk_Ndis,idx] = sort(Disk_Ndis);
    Disk_Nids = Disk_Nids(idx);
    Disk_theta = Disk_theta(idx);
    clear idx
    
    pick_ids = randsample(length(Disk_Nids),cut_num);
    Disk_Nids =  Disk_Nids(pick_ids);
    Disk_theta = Disk_theta(pick_ids);
    Disk_Ndis = Disk_Ndis(pick_ids);
    clear pick_ids
    
    [Disk_Ndis,idx] = sort(Disk_Ndis);
    Disk_Nids = Disk_Nids(idx);
    Disk_theta = Disk_theta(idx);
    clear idx
    
%     Disk_cors = N_X(Disk_Nids,:)';
%     figure; plot_mesh(X, F);hold on;
%     scatter3(X(Xid,1),X(Xid,2),X(Xid,3), 20, 'filled', 'r')
%     vis_ids_1 = find(Disk_theta>=0 & Disk_theta<pi);
%     vis_ids_2 = find(Disk_theta>=pi & Disk_theta<2*pi);
%     scatter3(Disk_cors(1,vis_ids_1),Disk_cors(2,vis_ids_1),Disk_cors(3,vis_ids_1), 2, 'filled', 'g');
%     scatter3(Disk_cors(1,vis_ids_2),Disk_cors(2,vis_ids_2),Disk_cors(3,vis_ids_2), 2, 'filled', 'y');
    
    %calculate Zernike bases in each disk
    r_k = Disk_Ndis/disk_R;
    theta_k = Disk_theta;
    Disk_bases = zernfun(Nn,Mm,r_k,theta_k,'norm');
    
    X_Disk_Nids(Xid,:)=Disk_Nids;
    X_Disk_Zerbases(Xid,:,:)=Disk_bases;
    
    %calculate input feature vec in each disk
    % 3d Height_map: aligned cordinate vector in disk (align Z-axis with
    % disk normal)
    Z_axis = [0;0;1];
    Disk_norm = center_normal;
    Disk_cors = N_X(Disk_Nids,:)';  
    vec_Disk_cors = Disk_cors - repmat(center_cors,1,length(Disk_cors));
    rot = vrrotvec(Disk_norm,Z_axis);
    mat = vrrotvec2mat(rot);
    Disk_FeaVec = mat*vec_Disk_cors;
    X_Disk_FeaVec(Xid,:,:)=Disk_FeaVec';
end

ZerPatch.resamplePids = X_Disk_Nids;
ZerPatch.bases = X_Disk_Zerbases;
ZerPatch.FeaVec = X_Disk_FeaVec;
end
