% process ADV data - PROFILE
% process data taken from multiple runs that form a single vertical profile
clear
close all

% plotting options
plot_raw_profs = 1;
plot_shift_profs = 1;
plot_phs_profs = 0;

if plot_raw_profs
    f_u_mean_raw = figure;
    f_restress_raw = figure;
end
if plot_shift_profs
    f_umean_shift = figure;
    f_uu_shift = figure;
end
f_Eu = figure;
f_Ew = figure;

N_bins_shift = 64;

% read metadata
warning off; 
addpath 'C:\Users\ljbak\My Drive\MATLAB\adv-velocity-matlab'
addpath 'C:\Users\ljbak\My Drive\MATLAB\fm-toolbox'
addpath 'H:\My Drive\MATLAB\fm-toolbox'

run_params = readtable('run_parameters.ods');
warning on; 

Np = 1; % profile number
runs = run_params.Run(run_params.profile==Np)'; 
[~,sort_idx] = sort(abs(run_params.ROI_z2_m(runs)));
runs = runs(sort_idx);

% preallocate
rng = [];
u = [];
v = [];
w1 = [];
w2 = [];
z = [];
eta_lr = [];

for n = runs
    fprintf('\nwind = %2.f Hz, current = %2.f Hz, %s\n', run_params.WindMotorFreq_Hz(n), ...
        run_params.CurrentPumpFreq_Hz(n), run_params.ADVDataFileName{n});
    
    %% read ADV data
    fn = [run_params.ADVDataFileName{n} '.mat'];
    load(fn);
   
    t = Data.Profiles_TimeStamp(1:run_params.ADVSampleCount(n));
    fs = 1/diff(t(1:2));
    rng_n = Data.Profiles_Range;
    rng = [rng;rng_n'];

    u_raw = Data.Profiles_VelX(1:run_params.ADVSampleCount(n),:);
    v_raw = Data.Profiles_VelY(1:run_params.ADVSampleCount(n),:);
    w1_raw = Data.Profiles_VelZ1(1:run_params.ADVSampleCount(n),:);
    w2_raw = Data.Profiles_VelZ2(1:run_params.ADVSampleCount(n),:);

    btm_t = Data.BottomCheck_TimeStamp;
    btm_z = Data.BottomCheck_BottomDistance;
        
    %% free surface
    % from bottom distance
    btm_z_mean = nanmean(btm_z);
    btm_z(abs(btm_z - btm_z_mean) > .1) = nan;
    btm_z_mean = nanmean(btm_z);
    A_rms = rms(btm_z - btm_z_mean,'omitnan');
    fprintf('mean surface distance = %1.3f cm, wave height = %1.3f cm\n',btm_z_mean*100, A_rms*100);
    
    % z coord w.r.t. water surface
    if run_params.ROI_z2_m(n)==0
        z = [z,repmat(rng_n,size(u_raw,1),1) - btm_z_mean]; 
        btm_z_ref = btm_z_mean;  % reference surface level for profile
    else
        z = [z,repmat(rng_n,size(u_raw,1),1) - btm_z_ref + run_params.ROI_z2_m(n)];
    end

%     % from amplitude
%     amp_thres = 60; % signal amplitude threshold for detecting surface
%     [surf_z, ~, ~, ~] = get_free_surface(amp, amp_thres, rng);
    
    % displacement relative to mean surface level (positive upward, negative downward)
%     eta_hr = surf_z - btm_z_mean;  % high-resolution eta
    
    eta_lr_n = interp1(btm_t,btm_z,t) - btm_z_mean; % low-resolution eta (interpolated to full resolution)

    tlim = 500;
%     figure; c = adv_prof_time_plot(t(1:tlim),rng,amp_mag(1:tlim,:));
%     title('Amp magnitude raw');
%     figure; c = adv_prof_time_plot(t(1:tlim),rng,amp_mag_sm(1:tlim,:));
%     title('Amp magnitude smoothed/filtered');
    
%     figure; plot(t(1:tlim),surf_z(1:tlim),'b-',btm_t(1:tlim/5),btm_z(1:tlim/5),'rx:'); 
%     xlabel('t [s]'); ylabel('free surface [m]'); xlim([t(1) t(tlim)])
%     legend('surface from amplitude','surface from bottom check')
    
    Ebar_norm = 1/(btm_t(end-10)-btm_t(10))*trapz(naninterp(btm_z(10:end-10) - btm_z_mean).^2);  % Ebar/(g*rho)
    H_rms = sqrt(8*Ebar_norm);
    
    if run_params.ROI_z2_m(n)==0
        eta_lr_nonan = eta_lr_n; eta_lr_nonan(isnan(eta_lr_nonan)) = []; % PDF/CDF of surface elevation
        [N_eta,edges_eta] = histcounts(eta_lr_nonan,100,'normalization','cdf');
        ctrs_eta = edges_eta(2:end) - diff(edges_eta)/2;
        [N_eta,ia] = unique(N_eta); ctrs_eta = ctrs_eta(ia);
        figure; histogram(abs(eta_lr_nonan),100); xlabel('\eta [m]'); ylabel('PDF')
        
        % dominant wave amplitude from 90th %ile of eta (since eta distribution
        % is cut off by adv settings - if had full eta distribution, would
        % choose 67th %ile to get significant wave height)
        A_dom = interp1(N_eta,ctrs_eta,.90); 
    end

    eta_lr = [eta_lr,repmat(eta_lr_n,1,size(u_raw,2))];
   

    %% velocity    
    
    % filter velocity: use PST (Parsheh N_bins_shift10) remove spikes
%     if exist([run_params.ADVDataFileName{n} '_despiked.mat'],'file')
%         load([run_params.ADVDataFileName{n} '_despiked.mat']);
%     else        
%         figs = 1;
%         for i = 1:size(u_raw,2)
%     %         if i==15; figs = 1; end
%             [UVW_new,Goodpoint] = mPST_ADVSpikeFilter(u_raw(:,i), v_raw(:,i), w1_raw(:,i), w2_raw(:,i), fs, figs);
%             u(:,i) = UVW_new(:,1);
%             v(:,i) = UVW_new(:,2);
%             w1(:,i) = UVW_new(:,3);
%             w2(:,i) = UVW_new(:,4);
%             clear UVW_new
%         end
%         [u_pdf, u_range] = pdf_var(u(:), 100, 0);
%         hold on; plot(u_range, u_pdf, 'r+');
%     
%         save([run_params.ADVDataFileName{n} '_despiked.mat'],'u','v','w1','w2');
% %     end
    

%     % correct velocities for ADV rotation
%     load(sprintf('R_run%i.mat',run_params.calibrationRun(n)))
% 
%     [UVW2_lab] = [u_raw(:),v_raw(:),w2_raw(:)]*R;
%     u_n = reshape(UVW2_lab(:,1), size(u_raw));
%     v_n = reshape(UVW2_lab(:,2), size(u_raw));
%     w2_n = reshape(UVW2_lab(:,3), size(u_raw));
%     [UVW1_lab] = [u_raw(:),v_raw(:),w1_raw(:)]*R;
%     w1_n = reshape(UVW1_lab(:,3), size(u_raw));

    u_n = u_raw;
    v_n = v_raw;
    w1_n = w1_raw;
    w2_n = w2_raw;

    % remove measurements above water surface
    echo_offset = 6e-3; % [m] distance below surface in which ADV velocities are contaminated by surface echo
    for i = 1:size(u_n,1)
        cutoff_rng = btm_z_ref + eta_lr(i) - echo_offset;
        u_n(i,rng_n >= cutoff_rng) = nan;
        v_n(i,rng_n >= cutoff_rng) = nan;
        w1_n(i,rng_n >= cutoff_rng) = nan;
        w2_n(i,rng_n >= cutoff_rng) = nan;
    end
    
    u = [u,u_n];
    v = [v,v_n];
    w1 = [w1,w1_n];
    w2 = [w2,w2_n];
end

%% raw mean profiles (function of z)
u_mean_raw = squeeze(nanmean(u,1));
v_mean_raw = squeeze(nanmean(v,1));
w1_mean_raw = squeeze(nanmean(w1,1));
w2_mean_raw = squeeze(nanmean(w2,1));

% plots
if plot_raw_profs
    % raw mean vel
    figure(f_u_mean_raw);
    h = adv_prof_subplots(z(1,:),[u_mean_raw', v_mean_raw', w1_mean_raw', w2_mean_raw'], ...
        'z [m]',{'\langle u\rangle [m/s]','\langle v\rangle [m/s]','\langle w_1\rangle [m/s]','\langle w_2\rangle [m/s]'});
%         for i = 1:length(h); set(h(i),'YLim',[.04 .07]); end
end

% raw velocity fluctuations
u_fluct_raw = u - repmat(u_mean_raw,size(u,1),1);
v_fluct_raw = v - repmat(v_mean_raw,size(v,1),1);
w1_fluct_raw = w1 - repmat(w1_mean_raw,size(w1,1),1);
w2_fluct_raw = w2 - repmat(w2_mean_raw,size(w2,1),1);

% raw Re stresses
uu_raw = u_fluct_raw.*u_fluct_raw;
vv_raw = v_fluct_raw.*v_fluct_raw;
w1w1_raw = w1_fluct_raw.*w1_fluct_raw;
w2w2_raw = w2_fluct_raw.*w2_fluct_raw;
uw1_raw = u_fluct_raw.*w1_fluct_raw;
uw2_raw = u_fluct_raw.*w2_fluct_raw;
w1w2_raw = w1_fluct_raw.*w2_fluct_raw;

uu_prof = nanmean(uu_raw,1);
vv_prof = nanmean(vv_raw,1);
w2w2_prof = nanmean(w2w2_raw,1);
uw2_prof = nanmean(uw2_raw,1);

if plot_raw_profs
    figure(f_restress_raw);
    h = adv_prof_subplots(z(1,:),[uu_prof', vv_prof', w2w2_prof', uw2_prof'], ...
        'z [m]',{'\langle u''u''\rangle [m/s]','\langle v''v''\rangle [m/s]','\langle w_2''w_2''\rangle [m/s]','\langle u''w_2''\rangle [m/s]'});
end

%% bricker and monismith wave-turb decomposition
%     if run_params.WindMotorFreq_Hz(n) > 0
%         uw1_wave = zeros(1,size(u,2));
%         for i = 1:size(u,2)
%             plot_on = 0;%i==15;
%             [U_wave,f_wave]=wave_turb_decomp(u(:,i),fs,plot_on);
%             [W1_wave,f_wave]=wave_turb_decomp(w1(:,i),fs,plot_on);
%             amp_U = abs(U_wave); 
%             amp_W1 = abs(W1_wave);
%             phs_U = atan2(imag(U_wave),real(U_wave));
%             phs_W = atan2(imag(W1_wave),real(W1_wave));
%             uw1_wave_f = abs(U_wave).*abs(W1_wave).*cos(phs_W - phs_U);
%             uw1_wave(i) = sum(uw1_wave_f);
%         end
% 
%         uw1_total = nanmean(u.*w1,1);
%         uw1_turb = uw1_total - uw1_wave;
%         figure; plot(uw1_total,rng,uw1_wave,rng,uw1_turb,rng)
%         legend('total','wave','turb'); xlabel('u''w'' [m/s]'); ylabel('z [m]'); title('decomposed stresses')
%         % switch to z-tilde coords?
% 
%         u_tau_bmdc = sqrt(abs(uw1_turb(end)));  % friction vel based on Bricker and Monismith decomposition
%         fprintf('u_tau_bm = %2.2f cm/s\n',u_tau_bmdc(end)*100);
%     end


%% velocity power spectra 
u_advec = u_mean_raw(1);%u_mean(end-5);  % advection velocity for Taylor's hyp/Doppler shift 

ss_idx = logical(rng>=.044 & rng<=.056); % index of sweet spot 

[E_u, f] = get_spectrum(u_fluct_raw(:,ss_idx), fs);
[E_w1, ~] = get_spectrum(w1_fluct_raw(:,ss_idx), fs);

figure(f_Eu);
loglog(f,E_u); hold on; 
if n == length(runs)
    [~,A_idx] = min(abs(f-1));
    A = E_u(A_idx);
    loglog(f(f>.1),A*f(f>.1).^(-5/3),'-k')
end
xlabel('f [s^{-1}]')
ylabel('E_u(f) [m^2/s^2/Hz]')

figure(f_Ew);
loglog(f,E_w1); hold on; 
if n == length(runs)
    [~,A_idx] = min(abs(f-1));
    A = E_w1(A_idx);
    loglog(f(f>.1),A*f(f>.1).^(-5/3),'-k')
end
xlabel('f [s^{-1}]')
ylabel('E_w(f) [m^2/s^2/Hz]')

[~,f_dom_idx] = max(E_w1.*(f>1)');  % maximum of vertical vel power spectrum
f_dom_o = f(f_dom_idx);   % observed (doppler-shifted) dominant wave freq 
f_dom = 1/(1/f_dom_o - -2*pi*u_advec/9.81);  % actual dominant wave freq (corrected for doppler shift)
cp_dom = 9.81/(2*pi*f_dom);   % dominant phase speed
k_dom = 2*pi*f_dom/cp_dom; %(cp_dom - u_advec);      % dominant wavenumber
L_dom = 2*pi/k_dom;  % dominant wavelength
eps_dom = A_rms*k_dom;  % dominant wave steepness
% eps_dom = run_params.WaveHeight_m(n)/2*k_dom;  % dominant wave steepness


%% velocity in free surface coords
% z positive upward, negative underwater, 0 at mean surface level  
% shifted z coord (z-tilde)
z_shift = z - eta_lr.*exp(k_dom*z);
%     z_shift = z - repmat(eta_lr,1,size(z,2)).*exp(k_dom*z);   

%     %% data points for which full wave orbital is available (based on 10th
%     % and 90th percentile surface elevation)
%     eta_c1 = interp1(N_eta, ctrs_eta, .9);
%     eta_c2 = interp1(N_eta, ctrs_eta, .1);
% 
%     z_shift_c1 = rng(end) - btm_z_mean - eta_c1;  
%     z_shift_c2 = rng(1) - btm_z_mean - eta_c2;
%     
%     if z_shift_c1 < z_shift_c2
%         z_shift_c2 = rng(1) - btm_z_mean;
%     end
% 
%     z_c1_c2_idx = z_shift(:) < z_shift_c1 & z_shift(:) > z_shift_c2;
z_c1_c2_idx = logical(~isnan(z_shift) & z_shift < 0);

%         figure; adv_prof_time_plot(t(1:tlim),rng,z_shift(1:tlim,:));
%         hold on; plot(t(1:tlim),(btm_z_mean+eta_lr(1:tlim))*1000,'k-')

% calculate phase using hilbert transform
u_fluct_bar = nanmean(u_fluct_raw,2);
%     for i = 1:size(u_fluct0_sm,2)
%         u_fluct0_sm(:,i) = smooth(u_fluct0_sm(:,i),5);
%     end
Y = hilbert(naninterp(u_fluct_bar));
phs = atan2(imag(Y),real(Y));
phs = repmat(phs,1,size(u,2));

% average over phase
N_phs = 12;
u_phs_mean = zeros(length(rng),N_phs);
for i = 1:length(rng)
    [u_phs_mean(i,:), phs_bins] = condition_vars(u(:,i),phs(:,i),N_phs,0);
end    
if plot_phs_profs
%     figure; pcolor(phs_bins*180/pi, rng, u_phs_mean'); shading flat; 
%     c = colorbar; c.Label.String = 'u [m/s]'; xlabel('\phi [deg]'); ylabel('z [m]')

    figure; 
    adv_phs_prof_plot(rng,u_phs_mean,phs_bins,'$z$ [m]','$\langleu(z, \phi)\rangle$ [m/s]');
end

% shifted mean velocity
[u_mean_shift,z_shift_prof,N] = condition_vars(u(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
v_mean_shift = condition_vars(v(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
w1_mean_shift = condition_vars(w1(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
w2_mean_shift = condition_vars(w2(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);

% shifted, phase-avg mean velocity
[u_mean_shift_phs,z_shift_prof,phs_bins,N] = condition_vars2(u(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
v_mean_shift_phs = condition_vars2(v(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
w1_mean_shift_phs = condition_vars2(w1(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
w2_mean_shift_phs = condition_vars2(w2(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);

if plot_phs_profs
    figure; 
    adv_phs_prof_plot(z_shift_prof,u_mean_shift_phs,phs_bins,'$\tilde{z}$ [m]','$\langle u(\tilde{z}, \phi)\rangle$ [m/s]');
end

if plot_shift_profs
    figure(f_umean_shift);
    adv_prof_subplots(z_shift_prof,[u_mean_shift, v_mean_shift, w1_mean_shift, w2_mean_shift], ...
        '$\tilde{z}$ [m]',{'\langle u\rangle [m/s]','\langle v\rangle [m/s]','\langle w_1\rangle [m/s]','\langle w_2\rangle [m/s]'});
end

% shifted velocity fluctuations
u_fluct_shift = u - interp1(z_shift_prof,u_mean_shift,z_shift);
v_fluct_shift = v - interp1(z_shift_prof,v_mean_shift,z_shift);
w1_fluct_shift = w1 - interp1(z_shift_prof,w1_mean_shift,z_shift);
w2_fluct_shift = w2 - interp1(z_shift_prof,w2_mean_shift,z_shift);

% instantaneous shifted Re stresses
uu_shift = u_fluct_shift.*u_fluct_shift;
vv_shift = v_fluct_shift.*v_fluct_shift;
vw2_shift = v_fluct_shift.*w1_fluct_shift;
w1w1_shift = w1_fluct_shift.*w1_fluct_shift;
w2w2_shift = w2_fluct_shift.*w2_fluct_shift;
uw1_shift = u_fluct_shift.*w1_fluct_shift;
uw2_shift = u_fluct_shift.*w2_fluct_shift;
w1w2_shift = w1_fluct_shift.*w2_fluct_shift;

% shifted Reynolds stress profiles 
uu_shift_prof = condition_vars(uu_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
vv_shift_prof = condition_vars(vv_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
w2w2_shift_prof = condition_vars(w2w2_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
uw2_shift_prof = condition_vars(uw2_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);
vw2_shift_prof = condition_vars(vw2_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),N_bins_shift,0);

if plot_shift_profs
    figure(f_uu_shift);
    adv_prof_subplots(z_shift_prof,[uu_shift_prof, vv_shift_prof, w2w2_shift_prof, vw2_shift_prof, uw2_shift_prof], ...
        '$\tilde{z}$ [m]',{'\langle u''u''\rangle [m/s]','\langle v''v''\rangle [m/s]','\langle w_2''w_2''\rangle [m/s]','\langle v''w_2''\rangle [m/s]','\langle u''w_2''\rangle [m/s]'});
end

%% shifted, phase-avg Reynolds stresses   
uu_shift_phs = condition_vars2(uu_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
vv_shift_phs = condition_vars2(vv_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
w2w2_shift_phs = condition_vars2(w2w2_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);
uw2_shift_phs = condition_vars2(uw2_shift(z_c1_c2_idx),z_shift(z_c1_c2_idx),phs(z_c1_c2_idx),[N_bins_shift,N_phs],[0 0]);

if plot_phs_profs
    figure; 
    adv_phs_prof_plot(z_shift_prof,uw2_shift_phs,phs_bins,'$\tilde{z}$ [m]','$\langle u''w_2''(\tilde{z}, \phi)\rangle$ [m/s]');
end


%     %% turbulence parameters
%     % stokes drift velocity 
%     if run_params.WindMotorFreq_Hz(n) > 0
%         U_S = 2*pi*f_dom*k_dom*(A_dom)^2;  % cp_dom*eps_dom; % 
%     else
%         U_S = 0;
%     end
%     fprintf('Stokes drift U_S = %2.2f cm/s\n',U_S*100);
%         
%     
%     %% shear velocity at free surface
%     % estimate from wind speed and C_D
%     vonkarman = 0.4;
%     rho_a = 1.2; % kg/m3
%     rho_w = 1000; % kg/m3
%     nu_w = 1e-6; % m2/s
% %     CD = 2.0e-3;  % drag coeff. 
%     z0 = 6.7e-4*(U10/cp_dom)^2.6*A_rms;    % *** FROM SAVELYEV ET AL 2020 ***
%     CD = vonkarman^2*(log(10/z0)).^(-2); 
%     z_anem = 0.8 - run_params.WaterDepth_m(n);
%     U10 = run_params.WindSpeed_m_s(n)*vonkarman/sqrt(CD)/(log(z_anem/10) + vonkarman/sqrt(CD));
%     u_tau_a_CD = sqrt(CD)*U10;
%     u_tau_CD = sqrt(rho_a/rho_w)*u_tau_a_CD;
%     fprintf('u*_a = %2.2f cm/s, u*_w = %2.2f cm/s\n',u_tau_a_CD*100,u_tau_CD*100);
%     
%    
%     % calculate u_tau from velocity gradient and shear stress
%     if run_params.WindMotorFreq_Hz(n) > 0
%         [~,dudz] = gradient(u_mean_shift_phs,diff(z_shift_prof(1:2)),diff(phs_bins(1:2)));
%         u_tau_momflux = sqrt(nu_w*abs(dudz) + abs(uw1_shift_phs));
%         figure; adv_phs_prof_plot(z_shift_prof,u_tau_momflux,phs_bins,'$\tilde{z}$ [m]','$u_\tau(\tilde{z}, \phi)\rangle$ [m/s]');
%     else
%         dudz = gradient(u_mean_raw,diff(z(1,1:2)));
%         u_tau_momflux = sqrt(nu_w*abs(dudz) + abs(uw1_prof));
%         figure; plot(u_tau_momflux,rng); xlabel('$u_\tau(\tilde{z})\rangle$ [m/s]','interp','latex'); ylabel('z [m]')
%     end
%     
%     
%     %% Langmuir number
%     La_t = sqrt(u_tau_CD/U_S);  % u_tau_bmdc  % max(u_tau_momflux(end,:))  % u_tau_CD
%     fprintf('La_t = %2.2f\n',La_t)

% end