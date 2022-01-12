function [R, U_lab, V_lab, W1_lab, varargout] = rotate_vels(U_adv, V_adv, W1_adv, W2_adv, varargin)
% [R, U_lab, V_lab, W1_lab, [W2_lab, R2]] = rotate_vels(U_adv, V_adv, W1_adv, W2_adv, [R])
% Takes ADV or ADV profiler velocity components in the instrument
% coordinate system and rotates them into laboratory coordinates, returning
% the rotated velocity components. Assumes laboratory x coordinate
% corresponds to the direction of maximum velocity; y and z are orthogonal
% to x.
%
% U_adv, etc: original velocities in ADV coordinates (NaNs in data not
% supported)
% U_lab, etc: rotated velocities in lab coordinates
% R: rotation matrix (optional)

% compute the SVD of velocities to determine the direction of maximum
% velocity
UVW_adv = [U_adv(:), V_adv(:), W1_adv(:)];
[L,S,R] = svd(UVW_adv, 'econ');
for i = 1:3
    if R(i,i) < 0
        R(:,i) = -R(:,i);
    end
end
UVW_lab = UVW_adv*R;

U_lab = reshape(UVW_lab(:,1), size(U_adv));
V_lab = reshape(UVW_lab(:,2), size(U_adv));
W1_lab = reshape(UVW_lab(:,3), size(U_adv));

% if this is an ADV profiler measurement, also compute the SVD using the
% W2 component 
if nargin > 3
    UVW2_adv = [U_adv(:), V_adv(:), W2_adv(:)];
    [L2,S2,R2] = svd(UVW2_adv, 'econ');
    for i = 1:3
        if R2(i,i) < 0
            R2(:,i) = -R2(:,i);
        end
    end
    UVW2_lab = UVW2_adv*R2;
    
    U_lab = reshape(UVW2_lab(:,1), size(U_adv));
    V_lab = reshape(UVW2_lab(:,2), size(U_adv));
    W2_lab = reshape(UVW2_lab(:,3), size(U_adv));
    varargout{1} = W2_lab;
    varargout{2} = R2;
    
%     % average the rotated U and V components computed with W1 and W2
%     U_lab = mean(cat(3, U_lab, U_lab2), 3);
%     V_lab = mean(cat(3, V_lab, V_lab2), 3);
end

end