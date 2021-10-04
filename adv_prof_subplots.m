function h = adv_prof_subplots(z,q,zlab,qlab)
% z: vertical coordinate (nz x 1 vector)
% q: quantities to plot as a function of z (nz x nq matrix)
% zlab: z axis label string (string)
% qlab: q axes label strings (nq x 1 cell containing strings)

nq = size(q,2);
h = zeros(nq,1);

for i = 1:nq
    h(i) = subplot(1,nq,i);
    adv_profile_plot(q(:,i),z); hold on 
    if i ~= 1
        ylabel(''); set(gca,'yTickLabel',[]);
    else
        ylabel(zlab,'interp','latex');
    end
    xlabel(qlab{i});
end
