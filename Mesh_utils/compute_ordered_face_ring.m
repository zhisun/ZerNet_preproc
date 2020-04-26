function vface = compute_ordered_face_ring(face_ring, vid)
% compute_ordered_face_ring  Reorder 1-ring faces in an adjacent order
%   INPUTS:
%       face_ring: 1-ring face neighborhood of vid. Matrix of 3 x (number of 1-ring faces).
%       vid: Vertex index. A scalar integer.
%   OUTPUT:
%       vface: Ordered 1-ring faces. Same dimension as face_ring
%
%   November 30, 2017 -- Code written by Stephen Baek (http://www.stephenbaek.com)

    vface = face_ring;
    
    % reorder vface such that the first vertex of vface 
    % is always the start vertex
    for i=1:size(vface,2)
        while vface(1,i)~=vid
            vface(:,i) = circshift(vface(:,i),1);
        end
    end
    
    % adjacency matrix indicating which faces of the vface are adjacent
    adj = spalloc(size(vface,2),size(vface,2), 2*size(vface,2));
    for i=1:size(vface,2)
        adj(i,find(vface(2,:)==vface(2,i) | vface(3,:)==vface(2,i) |...
                   vface(2,:)==vface(3,i) | vface(3,:)==vface(3,i)) ) = 1;
    end
    adj = adj - speye(size(adj));
    
    % check if vface is on boundary. If exists, a face on the boundary
    % comes to the first after reordering
    temp = vface([2,3],:);
    [temp, ia] = unique(sort(temp(:)));
    temp = find(diff(ia)==1);
    if ~isempty(temp)
        reorder = temp(1);
    else
        reorder = 1;
    end

    % reorder vface such that faces having adjacent index are neighbors
    idx = setdiff(1:size(vface,2),reorder);
    while ~isempty(idx)
        ic = setdiff(find(adj(reorder(end),:)),reorder);
        if ~isempty(ic)
            reorder = [reorder ic(1)];
        else
            reorder = reorder;
        end
        idx = setdiff(1:size(vface,2),reorder);
    end
    vface = vface(:,reorder);
    
    if vface(2,1) == vface(3,2)
        vface = vface(:,end:-1:1);
    end
end