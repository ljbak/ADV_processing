% calibration: find rotation matrix
clear
close all

% read metadata
warning off; 
addpath 'C:\Users\ljbak\My Drive\MATLAB\adv-velocity-matlab'
addpath 'C:\Users\ljbak\My Drive\MATLAB\fm-toolbox'
addpath 'H:\My Drive\MATLAB\fm-toolbox'

run_params = readtable('run_parameters.ods');
warning on; 

% run IDs corresponding to calibration cases
runs = [1]; % row vector

for n = runs
    fprintf('\nwind = %2.f Hz, current = %2.f Hz, %s\n', run_params.WindMotorFreq_Hz(n), ...
        run_params.CurrentPumpFreq_Hz(n), run_params.ADVDataFileName{n});
    
    % read ADV data
    fn = [run_params.ADVDataFileName{n} '.mat'];
    load(fn);
   
    rng = Data.Profiles_Range;
    ss_idx = logical(rng>=.044 & rng<=.056); % index of sweet spot 

    u_raw = Data.Profiles_VelX;
    v_raw = Data.Profiles_VelY;
    w1_raw = Data.Profiles_VelZ1;
    w2_raw = Data.Profiles_VelZ2;

    % get rotation matrix
    [~, ~, ~, ~, ~, R] = rotate_vels(u_raw(:,ss_idx), v_raw(:,ss_idx), w1_raw(:,ss_idx), w2_raw(:,ss_idx));
    [UVW2_lab] = [u_raw(:),v_raw(:),w2_raw(:)]*R;
    u_raw_r = reshape(UVW2_lab(:,1), size(u_raw));
    v_raw_r = reshape(UVW2_lab(:,2), size(u_raw));
    w2_raw_r = reshape(UVW2_lab(:,3), size(u_raw));

    % print and save
    disp(R)
    save(sprintf('R_run%i.mat',n),'R');

end