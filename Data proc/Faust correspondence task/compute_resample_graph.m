function resample_graph = compute_resample_graph(shape,sample_num)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
%
[resample_graph.N_X,resample_graph.I,resample_graph.B,resample_graph.r] = UniformSampling_on_mesh(X,F,sample_num,'Color','blue');
resample_graph.NumPids = sample_num;
end