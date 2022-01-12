% compare ADV velocity PDFs
clear
close all

nfigs = 3;
figs = zeros(nfigs,1);
for i = 1:nfigs
    figs(i) = figure;
end

% read metadata
warning off; 
addpath 'C:\Users\ljbak\My Drive\MATLAB\adv-velocity-matlab'
addpath 'C:\Users\ljbak\My Drive\MATLAB\fm-toolbox'
addpath 'H:\My Drive\MATLAB\fm-toolbox'

run_params = readtable('run_parameters.ods');
warning on; 

% enter nan to allow all values of that parameter
windsp = nan; % wind motor frequency
fetch = 6.8; % fetch
spanws = 0; % spanwise position
roi_z1 = -inf;%-0.12;
roi_z2 = inf;%-0.03;
prof_idx = (isnan(windsp)|(run_params.WindMotorFreq_Hz == windsp)) & ...
    (isnan(fetch)|(run_params.ROI_x_m == fetch)) & ...
    (isnan(spanws)|(run_params.ROI_y_m == spanws)) & ...
    run_params.WindMotorFreq_Hz > 0;
profs = unique(run_params.profile(prof_idx))'; % profile numbers

lg = 1; % include legend?

for n = profs
    load(sprintf('results_%02d.mat',n))

    z_idx = logical(z>roi_z1 & z<roi_z2);
    
    % u 
    u0 = u(:,z_idx);
    [u_pdf, u_range] = pdf_var(u0(:),100,0);
    figure(figs(1)); 
    plot(u_range,u_pdf,'.-','markersize',8);
    hold on;
    xlabel('u')

    % v
    u0 = v(:,z_idx);
    [u_pdf, u_range] = pdf_var(u0(:),100,0);
    figure(figs(2)); 
    plot(u_range,u_pdf,'.-','markersize',8);
    hold on;
    xlabel('v')
    
    % w2
    u0 = w2(:,z_idx);
    [u_pdf, u_range] = pdf_var(u0(:),100,0);
    figure(figs(3)); 
    plot(u_range,u_pdf,'.-','markersize',8);
    hold on;
    xlabel('w2')
end

for i = 1:nfigs
    figure(figs(i));
    goodplot([8 4])
end

%% legends
if lg
    lstr = cell(size(profs));
    for i = 1:length(lstr)
        idx = find(run_params.profile == profs(i));
        lstr{i} = sprintf('x=%1.1fm, \ny=%.2fm, \nU_w=%2.fm/s', ...
            run_params.ROI_x_m(idx(1)),run_params.ROI_y_m(idx(1)),run_params.WindSpeed_m_s(idx(1)));
    end

    for i = 1:nfigs
        figure(i);
        ac = flipud(get(gca,'Children'));
        legend(ac,lstr,'location','northeastoutside');
    end
end