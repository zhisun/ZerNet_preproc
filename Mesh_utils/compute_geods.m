function D = compute_geods(shape)
D = zeros(length(shape.X),length(shape.X));
shape.f_dns = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));
for vid = 1:length(shape.X)
    [~,D(:,vid)] = fast_marching(vid,shape,'vertex',0,1,shape.f_dns);
end

