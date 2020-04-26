function [Input_ZerCoeff, Recons_acc] = compute_Input_Zercoeff(shape,resample_graph,disk_R,n_rays,zernike_order,num_bases)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
D = shape.D_global;
D_N = resample_graph.D_N2N;
N = resample_graph.N_X;
I = resample_graph.I;
cloest_N2X = knnsearch(X,N);
%compute sampled surface points normal
[~,F_normal] = compute_normal(X,F);
N_normal = F_normal(:,I);
%
Nn = []; Mm = [];
for n = 0:zernike_order
    Nn = [Nn n*ones(1,n+1)];
    Mm = [Mm -n:2:n];
end

N_Disk_zernike_coeff = zeros(size(N,1),num_bases,size(N,2));
N_Disk_recons_acc = zeros(size(N,1),1);
%
vfring = compute_vertex_face_ring(F');

for Nid = 1:length(N)
    center_cors = N(Nid,:)';
    face_Vs = F(I(Nid),:);
    D_Nid = D_N(Nid,:);
    
    vid = cloest_N2X(Nid);
    D_vid = D(vid,:);
    facelist = vfring{vid};
    ring_Vs = F(facelist,:);
    ring_neibors_list = setdiff(unique(ring_Vs(:)),vid);
    ring_radius = mean(D_vid(ring_neibors_list));
    
    Disk_Nids = find(D_Nid < disk_R & D_Nid > 0.8*ring_radius);
    Disk_Nfaceids = I(Disk_Nids);
    center_faceid = I(Nid);
    Disk_Nids(Disk_Nfaceids==center_faceid)=[];
    
    Disk_Ndis = D_Nid(Disk_Nids);
    
    [directions, theta] = get_ray_directions(face_Vs, X', center_cors, n_rays);
    
    Disk_cors = N(Disk_Nids,:)';  
    vec_Disk_cors = Disk_cors - repmat(center_cors,1,length(Disk_cors)); 
    [~,Idxs]=pdist2(directions', vec_Disk_cors', 'cosine', 'Smallest', 1);
    Disk_theta = theta(Idxs);
    
    if (length(Disk_Nids) < num_bases)
        if (2*length(Disk_Nids)> num_bases)
            duplicate_num = num_bases-length(Disk_Nids);
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
            Disk_theta = [Disk_theta, THETA_add];
        else
            duplicate_ratio =  ceil((num_bases)/length(Disk_Nids));
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
    
    Disk_cors = N(Disk_Nids,:)';  
    [Disk_dis,tem_idx] = sort(Disk_Ndis);
    
    Disk_cors = Disk_cors(:,tem_idx);
    Disk_theta = Disk_theta(tem_idx);
    clear tem_idx
    
%     figure; plot_mesh(X, F);hold on;
%     scatter3(N(Nid,1),N(Nid,2),N(Nid,3), 20, 'filled', 'r')
%     vis_ids_1 = find(Disk_theta>=0 & Disk_theta<pi);
%     vis_ids_2 = find(Disk_theta>=pi & Disk_theta<2*pi);
%     scatter3(Disk_cors(1,vis_ids_1),Disk_cors(2,vis_ids_1),Disk_cors(3,vis_ids_1), 2, 'filled', 'g');
%     scatter3(Disk_cors(1,vis_ids_2),Disk_cors(2,vis_ids_2),Disk_cors(3,vis_ids_2), 2, 'filled', 'y');
    
    %calculate Zernike bases in each disk
    r_k = Disk_dis/disk_R;
    theta_k = Disk_theta;
    Disk_bases = zernfun(Nn,Mm,r_k,theta_k,'norm');
    
    %calculate input feature vec in each disk
    % 3d Height_map: aligned cordinate vector in disk (align Z-axis with
    % disk normal
    Z_axis = [0;0;1];
    Disk_norm = N_normal(:,Nid);
    vec_Disk_cors = Disk_cors - repmat(center_cors,1,length(Disk_cors));
    rot = vrrotvec(Disk_norm,Z_axis);
    mat = vrrotvec2mat(rot);
    Disk_FeaVec = mat*vec_Disk_cors;
    Disk_zernike_coeff = Disk_bases\Disk_FeaVec';
    
    recon_Disk_FeaVec = Disk_bases*Disk_zernike_coeff;
    recon_error =  abs((recon_Disk_FeaVec - Disk_FeaVec')./Disk_FeaVec');
    recons_acc = length(find(recon_error < 0.2))/length(recon_error(:));
        
    N_Disk_zernike_coeff(Nid,:,:) = Disk_zernike_coeff;
    N_Disk_recons_acc(Nid) = recons_acc;
end
Input_ZerCoeff = N_Disk_zernike_coeff;
Recons_acc = N_Disk_recons_acc;
end