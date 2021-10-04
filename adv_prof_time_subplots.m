function [c,h] = adv_prof_time_subplots(z,t,q,qlab)
% z: vertical coordinate (nz x 1 vector)
% t: time coordinate (nt x 1 vector)
% q: quantities to plot as a function of z and t (nz x nt x nq matrix)
% zlab: z axis label string (string)
% tlab: t axis label string (string)
% qlab: q axes label strings (nq x 1 cell containing strings)

nq = size(q,3);
h = zeros(nq,1);
c = zeros(nq,1);

for i = 1:nq
    h(i) = subplot(nq,1,i);
    c(i) = adv_prof_time_plot(t,z,q(:,:,i)); 
    title(qlab{i});
    if i ~= nq
        xlabel('');
    end
end