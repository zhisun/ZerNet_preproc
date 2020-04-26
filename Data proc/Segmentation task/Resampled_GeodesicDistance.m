function D_new = Resampled_GeodesicDistance(F,D,I,B)
D_new = zeros(length(I),length(I));
idx = F(I,:);
for i= 1:length(I)
    ids = F(I(i),:);
    D_i_toX = D(:,ids)*B(i,:)';
    temp_D = D_i_toX(idx);
    D_new(:,i)= dot(temp_D,B,2);
    D_new(i,i)= 0;
end
end