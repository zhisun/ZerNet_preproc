function [ray, theta] = get_ray_directions(face, V, center_cors, n_samples)
    
    v1 = V(:,face(1)) - center_cors;
    v2 = V(:,face(2)) - center_cors;
    v3 = V(:,face(3)) - center_cors;
    
    v1_normalized = v1./repmat(sqrt(sum(v1.^2)), 3,1);
    v2_normalized = v2./repmat(sqrt(sum(v2.^2)), 3,1);
    v3_normalized = v3./repmat(sqrt(sum(v3.^2)), 3,1);
    
    face_angle = [acos(dot(v1_normalized, v2_normalized)),acos(dot(v2_normalized, v3_normalized)),acos(dot(v3_normalized, v1_normalized))];
    angle_factor = 2*pi/sum(face_angle); % compensate angle defect (Gauss-Bonnet)
    face_angle = [face_angle; zeros(1, length(face_angle)); cumsum(face_angle*angle_factor)];
    face_angle(2,:) = [0, face_angle(3,1:end-1)];
        
    theta = ((0:n_samples-1)+0.5)/n_samples*2*pi;
%     theta = (0:n_samples-1)/n_samples*2*pi;
    
    idx = zeros(1, n_samples);
    mask = zeros(1, n_samples);
    for i=1:size(face_angle,2)
        flag = xor(theta < face_angle(3,i), mask);
        idx(find(flag)) = i;
        mask = mask | flag;
    end
    
    angle_offset = (theta - face_angle(2,idx));
    rotang = angle_offset ./ angle_factor;
    rotaxis = [cross(v1_normalized, v2_normalized),cross(v2_normalized, v3_normalized),cross(v3_normalized, v1_normalized)];
    rotaxis = rotaxis(:,idx);
    % Rodrigues rotation formula
    vs = [v1,v2,v3]; 
    ray = vs(:,idx).*repmat(cos(rotang),3,1) + cross(rotaxis, vs(:,idx)).*repmat(sin(rotang),3,1)...
         + rotaxis.*repmat(dot(rotaxis,vs(:,idx)),3,1).*repmat(1-cos(rotang),3,1);
end