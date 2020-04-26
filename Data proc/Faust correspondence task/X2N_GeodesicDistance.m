function D_X2N = X2N_GeodesicDistance(F,D,I,B)
D_X2N = zeros(length(D),length(I));
idx = F(I,:);
for Xid=1:length(D)
    D_Xid = D(:,Xid);
    D_tem = D_Xid(idx);
    D_Xid2N = dot(D_tem,B,2);
    D_X2N(Xid,:)= D_Xid2N';
end
end
