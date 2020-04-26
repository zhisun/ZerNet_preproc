function resample_graph = compute_resample_graph(shape,sample_num)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
J2m = shape.J2m;
%
[resample_graph.N_X,resample_graph.I,resample_graph.B,resample_graph.r] = UniformSampling_on_mesh(X,F,sample_num,'Color','blue');
resample_graph.NumPids = sample_num;
%compute interpolated stress on sampled surface points
tem_id = F(resample_graph.I,:);
tem_J2m = J2m(tem_id);
resample_graph.N_J2m = dot(tem_J2m,resample_graph.B,2);
end