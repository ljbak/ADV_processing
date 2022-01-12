% compare ADV velocity profiles
clearvars -except ct
close all 

nfigs = 6;
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

lg = 1; % include legend?

prof_idx = (isnan(windsp)|(run_params.WindMotorFreq_Hz == windsp)) & ...
    (isnan(fetch)|(run_params.ROI_x_m == fetch)) & ...
    (isnan(spanws)|(run_params.ROI_y_m == spanws)) & ...
    run_params.WindMotorFreq_Hz > 0;
profs = unique(run_params.profile(prof_idx))'; % profile numbers

for n = profs
    load(sprintf('results_%02d.mat',n))

    %% hi-res profiles
    % u mean abs
    figure(figs(1)); 
    h = adv_prof_subplots(z,[u_mean_abs'],..., v_mean_abs', w2_mean_abs'], ...
        '$z$ [m]',{'$\langle u\rangle$ [m/s]','$\langle v\rangle$ [m/s]','$\langle w_2\rangle$ [m/s]'}, lg);
    hold on;
    
    % re stress abs
    figure(figs(2));
    h = adv_prof_subplots(z(1,:),sqrt([uu_prof', vv_prof', w2w2_prof', uw2_prof']), ...
        '$z$ [m]',{'$\sqrt{\langle u''u''\rangle}$ [m/s]','$\sqrt{\langle v''v''\rangle}$ [m/s]','$\sqrt{\langle w_2''w_2''\rangle}$ [m/s]', ...
        '$\sqrt{\langle u''w_2''\rangle}$ [m/s]'}, lg);
    hold on;
    
    % u mean shifted
    figure(figs(3));
    adv_prof_subplots(z_shift_prof,[u_mean_shift],..., v_mean_shift, w2_mean_shift], ...
        '$\tilde{z}$ [m]',{'$\langle u\rangle$ [m/s]','$\langle v\rangle$ [m/s]','$\langle w_2\rangle$ [m/s]'}, lg);
    hold on;

    % re stress shifted
    figure(figs(4));
    adv_prof_subplots(z_shift_prof,sqrt([uu_shift_prof, vv_shift_prof, w2w2_shift_prof, uw2_shift_prof]), ...
        '$\tilde{z}$ [m]',{'$\sqrt{\langle u''u''\rangle}$ [m/s]','$\sqrt{\langle v''v''\rangle}$ [m/s]','$\sqrt{\langle w_2''w_2''\rangle}$ [m/s]', ...
        '$\sqrt{\langle u''w_2''\rangle}$ [m/s]'}, lg);
    hold on;

    %% lo-res profiles
    z_lr = 0.015:-.03*4/5:min(z);
    
    % u mean abs
    figure(figs(5)); 
    h = adv_prof_subplots(z_lr,[interp1(z,u_mean_abs,z_lr)', interp1(z,v_mean_abs,z_lr)', interp1(z,w2_mean_abs,z_lr)'], ...
        '$z$ [m]',{'$\langle u\rangle$ [m/s]','$\langle v\rangle$ [m/s]','$\langle w_2\rangle$ [m/s]'}, lg);
    hold on;
    
    % re stress abs
    figure(figs(6));
    h = adv_prof_subplots(z_lr,[interp1(z,uu_prof,z_lr)', interp1(z,vv_prof,z_lr)', interp1(z,w2w2_prof,z_lr)', interp1(z,vw2_prof,z_lr)', interp1(z,uw2_prof,z_lr)'], ...
        '$z$ [m]',{'$\langle u''u''\rangle$ [m/s]','$\langle v''v''\rangle$ [m/s]','$\langle w_2''w_2''\rangle$ [m/s]', ...
        '$\langle v''w_2''\rangle$ [m/s]','$\langle u''w_2''\rangle$ [m/s]'}, lg);
    hold on;
    
    %% Log layer profile
%     ut = 0.01; %0.0148; %0.0181; % friction velocity estimate
%     nu = 1e-6;
%     del_nu = nu/ut;    % viscous length 
%     k = .41; B = 5.5; %5.2;  % assume canonical log law coefficients
%     rhs = 1/k*log(-z_shift_prof/del_nu)+B;
%     % rhs2 = y/del_nu;
% 
%     u_ll = -u_mean_shift + u_mean_shift(end);
%     
%     figure; semilogx(-z_shift_prof/del_nu,u_ll/ut,'k.',-z_shift_prof/del_nu,rhs,'r--'); hold on
%     % plot(y(3:5)/del_nu,rhs2(3:5),'b--')
%     xlabel('$-\tilde{z}^+$','interpreter','latex'); ylabel('u^+'); title('Mean velocity profile'); grid on; %goodplot([7 6]);
%     fprintf(['u_tau = ' num2str(ut) ' m/s\n'])
%     keyboard

end

% legends
if lg
    lstr = cell(size(profs));
    for i = 1:length(lstr)
        idx = find(run_params.profile == profs(i));
%         lstr{i} = sprintf('$x=%1.1f$m, \n$y=%.2f$m, \n$U_w=%2.f$m/s', ...
%             run_params.ROI_x_m(idx(1)),run_params.ROI_y_m(idx(1)),run_params.WindSpeed_m_s(idx(1)));
        lstr{i} = sprintf('$U_w=%2.f$m/s', run_params.WindSpeed_m_s(idx(1)));
    end

    for i = 1:nfigs
        figure(i);
        fc = get(gcf,'Children');
        ac = flipud(get(fc(2),'Children'));
        legend(ac,lstr,'position',[.8 .4 .1 .1]);
    end
end

for i = 1:nfigs
    figure(figs(i));
    goodplot([8 4])
end