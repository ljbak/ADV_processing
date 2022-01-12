% make z-t velocity plots
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
windsp = 60; % wind motor frequency
fetch = 6.8; % fetch
spanws = 0; % spanwise location
roi_z2 = -.06; % depth of top of sampling volume

runs = run_params.Run( (isnan(windsp)|(run_params.WindMotorFreq_Hz == windsp)) & ...
    (isnan(fetch)|(run_params.ROI_x_m == fetch)) & ...
    (isnan(spanws)|(run_params.ROI_y_m == spanws)) & ...
    (isnan(roi_z2)|(run_params.ROI_z2_m == roi_z2)) )';

for n = runs
    % load and preprocess data
    fn = [run_params.ADVDataFileName{n} '.mat'];
    load(fn);
    
    t = Data.Profiles_TimeStamp;
    fs = 1/diff(t(1:2));
    rng = Data.Profiles_Range;
    
    u_raw = Data.Profiles_VelX;
    v_raw = Data.Profiles_VelY;
    w1_raw = Data.Profiles_VelZ1;
    w2_raw = Data.Profiles_VelZ2;

    btm_t = Data.BottomCheck_TimeStamp;
    btm_z = Data.BottomCheck_BottomDistance;
    
    % z from bottom distance
    btm_z_mean = nanmean(btm_z);
    btm_z(abs(btm_z - btm_z_mean) > .1) = nan;
    btm_z_mean = nanmean(btm_z);
    A_rms = rms(btm_z - btm_z_mean,'omitnan');
    fprintf('mean surface distance = %1.3f cm, wave height = %1.3f cm\n',btm_z_mean*100, A_rms*100);
    
    % z coord w.r.t. water surface
    z = repmat(rng,size(u_raw,1),1) - btm_z_mean;  
    eta_lr = interp1(btm_t,btm_z,t) - btm_z_mean; % low-resolution eta (interpolated to full resolution)
    
    if exist(sprintf('R_%03d.mat',run_params.calibrationRun(n)),'file')
        load(sprintf('R_%03d.mat',run_params.calibrationRun(n)))
        [UVW2_lab] = [u_raw(:),v_raw(:),w2_raw(:)]*R;
        u = reshape(UVW2_lab(:,1), size(u_raw));
        v = reshape(UVW2_lab(:,2), size(u_raw));
        w2 = reshape(UVW2_lab(:,3), size(u_raw));
    else
        u = u_raw;
        v = v_raw;
        w2 = w2_raw;
    end
    if run_params.ROI_y_m(n) >= 0
        v = -v;
        u = -u;
    end
    echo_offset = 6e-3; % [m] distance below surface in which ADV velocities are contaminated by surface echo
    for i = 1:size(u,1)
        cutoff_rng = btm_z_mean + eta_lr(i) - echo_offset;
        u(i,rng >= cutoff_rng) = nan;
        v(i,rng >= cutoff_rng) = nan;
        w2(i,rng >= cutoff_rng) = nan;
    end
    %%
    figure;
    tstr = sprintf('$x=%1.1f$m, $y=%.2f$m, $U_w=%2.f$m/s, $%1.2f<z<%1.2f$m', ...
            run_params.ROI_x_m(n),run_params.ROI_y_m(n),run_params.WindSpeed_m_s(n),...
            run_params.ROI_z1_m(n),run_params.ROI_z2_m(n));
    adv_prof_time_subplots(z(1,:),t,cat(3,u,v,w2),{'u [m/s]','v [m/s]','w_2 [m/s]'},tstr);    

    goodplot([8 6])
end