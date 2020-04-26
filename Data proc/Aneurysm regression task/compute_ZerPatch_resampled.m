function ZerPatch = compute_ZerPatch_resampled(shape, resample_graph, disk_R, n_rays, cut_num, zernike_order, num_bases, min_r)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
D_N2N = Resampled_GeodesicDistance(F,shape.D_global,resample_graph.I,resample_graph.B);
D_N = D_N2N;
N = resample_graph.N_X;
I = resample_graph.I;
%
vfring = compute_vertex_face_ring(F');
shape.idxs = vfring;
shape.f_dns = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));

%compute sampled surface points normal
[~,F_normal] = compute_normal(X,F);
N_normal = F_normal(:,I);
%
Nn = []; Mm = [];
for n = 0:zernike_order
    Nn = [Nn n*ones(1,n+1)];
    Mm = [Mm -n:2:n];
end
%
N_Disk_Nids = zeros(length(N),cut_num);
N_Disk_Zerbases = zeros(length(N),cut_num, num_bases);
N_Disk_FeaVec = zeros(length(N),cut_num, 3);
%

for Nid = 1:length(N)
    
    center_cors = N(Nid,:)';
    D_Nid = D_N(Nid,:);
    
    %
    center_normal = N_normal(:,Nid);
    [~, tem] = mink(D_Nid,2);
    nearest_Nid = tem(2);
    
    ref_dir = N(nearest_Nid,:)'-center_cors;
    projected_ref_dir = ref_dir - center_normal*dot(ref_dir,center_normal);
    projected_ref_dir = projected_ref_dir/norm(projected_ref_dir);
    
    same_face_Nids = find(I==I(Nid));
    Dis = vecnorm(N(same_face_Nids,:)'-repmat(center_cors,1,length(same_face_Nids)));
    D_Nid(same_face_Nids) = Dis;
    
    lower_bound = max(vecnorm(N(same_face_Nids,:)'-repmat(center_cors,1,length(same_face_Nids))));
    lower_bound = min(lower_bound, min_r);
    
    Neibor_Nids = find(D_Nid < disk_R & D_Nid > lower_bound);
    Neibor_Ndis = D_Nid(Neibor_Nids);
    Neibor_dirs = N(Neibor_Nids,:)'-repmat(center_cors,1,length(Neibor_Nids));
    
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
    
%     Disk_cors = N(Disk_Nids,:)';
%     figure; plot_mesh(X, F);hold on;
%     scatter3(N(Nid,1),N(Nid,2),N(Nid,3), 20, 'filled', 'r')
%     vis_ids_1 = find(Disk_theta>=0 & Disk_theta<pi);
%     vis_ids_2 = find(Disk_theta>=pi & Disk_theta<2*pi);
%     scatter3(Disk_cors(1,vis_ids_1),Disk_cors(2,vis_ids_1),Disk_cors(3,vis_ids_1), 2, 'filled', 'g');
%     scatter3(Disk_cors(1,vis_ids_2),Disk_cors(2,vis_ids_2),Disk_cors(3,vis_ids_2), 2, 'filled', 'y');
    
    %calculate Zernike bases in each disk
    r_k = Disk_Ndis/disk_R;
    theta_k = Disk_theta;
    Disk_bases = zernfun(Nn,Mm,r_k,theta_k,'norm');
    
    N_Disk_Nids(Nid,:)=Disk_Nids;
    N_Disk_Zerbases(Nid,:,:)=Disk_bases;
    
    %calculate input feature vec in each disk
    % 3d Height_map: aligned cordinate vector in disk (align Z-axis with
    % disk normal)
    
    Z_axis = [0;0;1];
    Disk_norm = center_normal;
    Disk_cors = N(Disk_Nids,:)';  
    vec_Disk_cors = Disk_cors - repmat(center_cors,1,length(Disk_cors));
    rot = vrrotvec(Disk_norm,Z_axis);
    mat = vrrotvec2mat(rot);
    Disk_FeaVec = mat*vec_Disk_cors;
    N_Disk_FeaVec(Nid,:,:)=Disk_FeaVec';
end

ZerPatch.resamplePids = N_Disk_Nids;
ZerPatch.bases = N_Disk_Zerbases;
ZerPatch.FeaVec = N_Disk_FeaVec;
end
