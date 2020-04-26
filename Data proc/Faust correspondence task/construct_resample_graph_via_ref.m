function resample_graph = construct_resample_graph_via_ref(shape,ref_resample_graph)
X = [shape.X,shape.Y,shape.Z];
F = shape.TRIV;
%
resample_graph.I = ref_resample_graph.I;
resample_graph.B = ref_resample_graph.B;
resample_graph.NumPids = ref_resample_graph.NumPids;

%compute sampled surface points via interpolation weights on original mesh vertices
tem_id = F(resample_graph.I,:);
N_X = zeros(length(resample_graph.I),3);
for i =1:3
    V1 = X(tem_id(:,1),:);
    V2 = X(tem_id(:,2),:);
    V3 = X(tem_id(:,3),:);
    N_X(:,i) = dot([V1(:,i),V2(:,i),V3(:,i)],resample_graph.B,2);
end
resample_graph.N_X = N_X;
end