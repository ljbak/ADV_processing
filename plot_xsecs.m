% plot a velocity cross-section
clear
close all

% read metadata
warning off; 
addpath 'C:\Users\ljbak\My Drive\MATLAB\adv-velocity-matlab'
addpath 'C:\Users\ljbak\My Drive\MATLAB\fm-toolbox'
addpath 'H:\My Drive\MATLAB\fm-toolbox'

run_params = readtable('run_parameters.ods');
warning on; 

% enter nan to allow all values of that parameter
windsp = 45; % wind motor frequency
fetch = 6.8; % fetch

prof_idx = (isnan(windsp)|(run_params.WindMotorFreq_Hz == windsp)) & ...
    (isnan(fetch)|(run_params.ROI_x_m == fetch)) & ...
    run_params.WindMotorFreq_Hz > 0;
profs = unique(run_params.profile(prof_idx))'; % profile numbers

%% 3d
figure;

for n = profs
    load(sprintf('results_%02d.mat',n))

    spanws = run_params.ROI_y_m(run_params.profile==n);
    [zz, yy, xx] = meshgrid(z_shift_prof,spanws(1),fetch);
    quiver3(xx,yy,zz,u_mean_shift',v_mean_shift',w2_mean_shift');
    hold on
    
end

axis equal
axis vis3d
ylim([-.45 .45])
zlim([-.6 0.1])

%% 2d
figure;

for n = profs
    load(sprintf('results_%02d.mat',n))

    spanws = run_params.ROI_y_m(run_params.profile==n);
    [zz, yy] = meshgrid(z_shift_prof,spanws(1));
    quiver(yy,zz,v_mean_shift',w2_mean_shift');
    hold on
    
end

axis equal
axis([-.45 .45 -.6 0.1])