function [split_num,mean_disk_num] = compute_disk_num(resample_graph,disk_R)
X_Disk_Nids_num = zeros(size(resample_graph.D_X2N,1),1);
for Xid = 1:length(X_Disk_Nids_num)
    D_Xid2N = resample_graph.D_X2N(Xid,:);
    Disk_Nids = find(D_Xid2N < disk_R & D_Xid2N > 0);
    X_Disk_Nids_num(Xid)= length(Disk_Nids);
end
split_num = ceil(mean(X_Disk_Nids_num) - std(X_Disk_Nids_num));
mean_disk_num = floor(mean(X_Disk_Nids_num));
end