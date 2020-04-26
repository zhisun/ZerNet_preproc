function [ray, theta] = get_sample_directions(vertex, face, vfring, vid, n_samples)
% get_sample_directions  Sample (n_samples) rays eminating from the vertex (vid)
%   INPUTS:
%       vertex: A 3 by (number of vertices) matrix containing the vertex coordinates.
%       face: A 3 by (number of faces) matrix containing the triangle indices
%       vfring: Vertex face ring. Cell array of dimension (number of vertices).
%               Each cell contains a face index array neighboring to the corresponding vertex.
%       vid: Vertex index. A scalar integer.
%       n_samples: Number of rays to be sampled.
%   OUTPUT:
%       ray: A 3 by n_samples matrix whose column vectors are the sampled rays.
%            Rays are sampled in counter-clock-wise order.
%
%   November 30, 2017 -- Code written by Stephen Baek (http://www.stephenbaek.com)
    vface = compute_ordered_face_ring(face(:, vfring{vid}), vid);
    
    v1 = vertex(:,vface(2,:)) - vertex(:,vface(1,:));
    v2 = vertex(:,vface(3,:)) - vertex(:,vface(1,:));
    
    v1_normalized = v1./repmat(sqrt(sum(v1.^2)), 3,1);
    v2_normalized = v2./repmat(sqrt(sum(v2.^2)), 3,1);
    
    face_angle = acos(dot(v1_normalized, v2_normalized));
    angle_factor = 2*pi/sum(face_angle); % compensate angle defect (Gauss-Bonnet)
    face_angle = [face_angle; zeros(1, length(face_angle)); cumsum(face_angle*angle_factor)];
    face_angle(2,:) = [0, face_angle(3,1:end-1)];
        
    theta = ((0:n_samples-1)+0.5)/n_samples*2*pi;
    
    idx = zeros(1, n_samples);
    mask = zeros(1, n_samples);
    for i=1:size(face_angle,2)
        flag = xor(theta < face_angle(3,i), mask);
        idx(find(flag)) = i;
        mask = mask | flag;
    end
    
    angle_offset = (theta - face_angle(2,idx));
    rotang = angle_offset ./ angle_factor;
    rotaxis = cross(v1_normalized, v2_normalized);
    rotaxis = rotaxis(:,idx);
    
    % Rodrigues rotation formula
    ray = v1(:,idx).*repmat(cos(rotang),3,1) + cross(rotaxis, v1(:,idx)).*repmat(sin(rotang),3,1)...
         + rotaxis.*repmat(dot(rotaxis,v1(:,idx)),3,1).*repmat(1-cos(rotang),3,1);
end