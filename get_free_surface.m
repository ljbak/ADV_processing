function [surf_z, surf_idx, amp_mag, amp_mag_sm] = get_free_surface(amp, amp_thres, rng)
% get free surface distance surf_z and corresponding index surf_idx from
% ADV profiler signal amplitude amp and range rng, considering only regions
% where the amplitude is below a value of amp_thres. rng(surf_idx(i))
% corresponds to eta(i) (procedure from Chemin, Luneau, & Caulliez 2020)

% magnitude of signal amplitude
amp_mag = (amp(:,:,1).^2+amp(:,:,2).^2+amp(:,:,3).^2+amp(:,:,4).^2).^(1/2);

% preallocate
amp_mag_sm = amp_mag;
surf_z = zeros(size(amp_mag,1),1);
surf_idx = zeros(size(amp_mag,1),1);

for i = 1:size(amp_mag,1)
    % smooth ADV signal amplitude
    amp_mag_sm(i,:) = smooth(amp_mag(i,:),3);
    
    % look for next eta within search_dist of previous eta
    search_dist = .005; 
    if i > 1
        if ~isnan(surf_z(i-1))
%                 amp_mag_sm(i,rng < (eta(i-1) - search_dist) | rng > (eta(i-1) + search_dist)) = nan;
%             else
%                 amp_mag_sm(i,rng < (eta(i-1) - search_dist) | rng > (eta(i-1) + search_dist)) = nan;
        end
    else
        amp_mag_sm(i,1:floor(2*end/3)) = nan;
    end
    
    amp_mag_fit = amp_mag_sm;  % preserve all values for parabolic fit later
    amp_mag_sm(i,amp_mag_sm(i,:) > amp_thres) = nan;  % nan values that are outside search distance
    
    if any(~isnan(amp_mag_sm(i,:)))
        [val(i), I] = min(amp_mag_sm(i,:));  % find amplitude minimum
        surf_idx(i) = I;  
        
        % 3-pt parabolic fit for eta (if there are data points on either side of eta)
        surf_z(i) = rng(I);
        if I < length(rng) && I > 1
            coeffs = [rng(I-1:I+1)'.^2 rng(I-1:I+1)' [1;1;1]]\amp_mag_fit(i,I-1:I+1)';
            surf_z(i) = -coeffs(2)/(2*coeffs(1));
        end
    else
        surf_idx(i) = nan;
        surf_z(i) = nan;
    end
end

end