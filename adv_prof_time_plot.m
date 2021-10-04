function c = adv_prof_time_plot(t,rng,q)
% plot variable q along the ADV profile range over time

pcolor(t,rng*1000,q'); 
shading flat; 
colormap jet
xlabel('time [s]')
ylabel('z [mm]')
c = colorbar;