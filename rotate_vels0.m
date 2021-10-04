function [U_lab, V_lab, W1_lab, varargout] = rotate_vels0(U_adv, V_adv, W1_adv, W2_adv)
% Takes ADV or ADV profiler velocity components in the instrument
% coordinate system and rotates them into laboratory coordinates, returning
% the rotated velocity components. Assumes laboratory x coordinate
% corresponds to the direction of maximum velocity; y and z are orthogonal
% to x.
%
% U_adv, etc: original velocities in ADV coordinates (NaNs in data not
% supported)
% U_lab, etc: rotated velocities in lab coordinates

% compute the SVD of velocities to determine the direction of maximum
% velocity
[L,S,R] = svd([U_adv(:), V_adv(:), W1_adv(:)], 'econ');
% L = L.*sign(R); % FIX SIGNS - generalize from current-only case
UVW = L*S*sign(R);

U_lab = reshape(UVW(:,1), size(U_adv));
V_lab = reshape(UVW(:,2), size(U_adv));
W1_lab = reshape(UVW(:,3), size(U_adv));

% if this is an ADV profiler measurement, also compute the SVD using the
% W2 component 
if nargin > 3
    
    [L2,S2,R2] = svd([U_adv(:), V_adv(:), W2_adv(:)], 'econ');
%     L2 = L2.*sign(R2);
    UVW2 = L2*S2;
    
    U_lab2 = reshape(UVW2(:,1), size(U_adv));
    V_lab2 = reshape(UVW2(:,2), size(U_adv));
    W2_lab = reshape(UVW2(:,3), size(U_adv));
    varargout{1} = W2_lab;
    
    % average the rotated U and V components computed with W1 and W2
    U_lab = mean(cat(3, U_lab, U_lab2), 3);
    V_lab = mean(cat(3, V_lab, V_lab2), 3);
end

end