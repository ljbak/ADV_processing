function adv_phs_prof_plot(z,q,phs_bins,zlab,qlab)
% z: vertical coordinate (nz x 1 vector)
% phs_bins: phase bins (N_phs x 1 vector)
% q: phase-binned quantity to plot as a function of z (nz x N_phs matrix)
% zlab: z axis label string (string)
% qlab: q axis label string (string)

N_phs = length(phs_bins);

c = parula(N_phs);
lg = cell(N_phs,1);
for i = 1:N_phs
    plot(q(:,i),z','color',c(i,:),'linewidth',1.25); hold on; 
    lg{i} = ['\phi=' num2str(phs_bins(i)*180/pi,3) '\circ'];
end
ylabel(zlab,'interp','latex'); xlabel(qlab,'interp','latex');
legend(lg,'location','northeastoutside');
