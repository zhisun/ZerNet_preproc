function N_Disk_Nids_num = compute_disk_num(resample_graph,disk_R)
% X_Disk_Nids_num = zeros(size(resample_graph.D_X2N,1),1);
% for Xid = 1:length(X_Disk_Nids_num)
%     D_Xid2N = resample_graph.D_X2N(Xid,:);
%     Disk_Nids = find(D_Xid2N < disk_R & D_Xid2N > 0);
%     X_Disk_Nids_num(Xid)= length(Disk_Nids);
% end
% split_num = ceil(mean(X_Disk_Nids_num) - std(X_Disk_Nids_num));
% mean_disk_num = floor(mean(X_Disk_Nids_num));

N_Disk_Nids_num = zeros(size(resample_graph.D_N2N,1),1);
for Nid = 1:length(N_Disk_Nids_num)
    D_Nid2N = resample_graph.D_N2N(Nid,:);
    Disk_Nids = find(D_Nid2N < disk_R & D_Nid2N > 0);
    
    Disk_Nfaceids = resample_graph.I(Disk_Nids);
    center_faceid = resample_graph.I(Nid);
    Disk_Nids(Disk_Nfaceids==center_faceid)=[];
    
    N_Disk_Nids_num(Nid)= length(Disk_Nids);
end
end