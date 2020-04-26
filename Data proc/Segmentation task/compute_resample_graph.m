function resample_graph = compute_resample_graph(shape,sample_num)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
% D = shape.D_global;
label_X = shape.label_V;
label_F = shape.label_F;
%
[resample_graph.N_X,resample_graph.I,resample_graph.B,resample_graph.r] = UniformSampling_on_mesh(X,F,sample_num,'Color','blue');
resample_graph.NumPids = sample_num;
% resample_graph.D_X2N = X2N_GeodesicDistance(F,D,resample_graph.I,resample_graph.B);
% resample_graph.D_N2N = Resampled_GeodesicDistance(F,D,resample_graph.I,resample_graph.B);

%compute labels on sampled surface points
resample_graph.N_label = label_F(resample_graph.I);
end