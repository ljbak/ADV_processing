% automated data-taking function
% https://github.com/jason-holloway/velmex-stage-control
%
% before running this program: 
%
% 1. HOME BISLIDE: manually adjust horizontal bislide position to center
% spanwise location and vertical axis so that the transducer is 7 cm below
% still surface level (i.e. top of sampling volume is at the still surface
% level). this is home (origin). note: positive displacements move carriage
% away from the motor. for y<0 cases: home is when transducer is at
% y = -20cm.
%
% 2. SET UP ADV SOFTWARE: Set up ADV in Nortek software and align top left
% corner of software wind with top left corner of screen. This allows the
% record button to be in a predictable location.
%
% 3. SET RUN_PARAMETERS FILEPATH, WINDSPEED, SPANWISE SIGN, AND FETCH:
% update the filepath of the run_parameters.ods spreadsheet, windspeed,
% spanwise sign, and fetch for the current set of runs (below)

addpath('C:\Users\Naeem\Documents\MATLAB\velmex-stage-control-master')

%% motor set-up
stage = VMX(3,1000,[],1,[],2)
stage.setCurrentPositionAsHome('y');
stage.setCurrentPositionAsHome('z');

%% inputs for each set of runs
windspeed = 45;  % speed of wind motor for this set of runs (Hz)
spanwise_sign = 1; % +1 for side of channel away from the horz motor (including centerline); -1 for side of channel near the horz motor
fetch = 6.8;  % fetch for this set of runs (m)

% cursor set-up and execute experiment
% cursor control and adv record button
hroot = groot;          % cursor object
% disp(hroot.PointerLocation);  % move cursor over record button and run this line to get button location
ADV_start_button_loc = [42 1092];  % position on screen of ADV start button (px)
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;

buffer_time = 60;       % time to wait between end of recording and start of bislide move (s)

% ADV run parameters 
warning off; 
run_params = readtable('G:\My Drive\MP in OSBL\luci adv data\run 211209\run_parameters.ods');
warning on;

% select the set of runs corresponding to the current wind speed and fetch
if spanwise_sign == 1
    runs = find(run_params.WindMotorFreq_Hz==windspeed & run_params.ROI_x_m==fetch & run_params.ROI_y_m>=0);
    y_offset = 0;
else
    runs = find(run_params.WindMotorFreq_Hz==windspeed & run_params.ROI_x_m==fetch & run_params.ROI_y_m<0);
    y_offset = 200;
end

% % starting after fault
% runs = runs(runs>3);

N_runs = length(runs);

try
    % move the bislide motors and start ADV data collection
    for n = runs'
        fprintf('\nRun %i (%i/%i)\n',run_params.Run(n),n-runs(1)+1,N_runs)

        % move bislide to next sampling position & wait for completion
        fprintf('moving bislide...')
        stage.moveMotorAbsolute('y',run_params.ROI_y_m(n)*1000 + y_offset); % y axis move (mm)
        stage.moveMotorAbsolute('z',-run_params.ROI_z2_m(n)*1000); % z axis move (mm) (switch sign)
        fprintf('complete\n')
        
        % move cursor to ADV record button
        pause(1)
        set(hroot,'PointerLocation',ADV_start_button_loc);
        pause(1)

        % click ADV record button in Nortex software
        mouse.mousePress(InputEvent.BUTTON1_MASK);    % left click press
        mouse.mouseRelease(InputEvent.BUTTON1_MASK);  % left click release

        % wait for recording to finish
        fprintf('recording ADV data...')
        pause(run_params.ADVDuration_min(n)*60 + buffer_time);
        fprintf('complete\n')

    end
catch
    fprintf('\nerror on Run %i\n',n)
    return
end

%% rehome and take bislide offline
fprintf('\nRehoming\n')
stage.moveMotorAbsolute('y',0);
stage.moveMotorAbsolute('z',0);

str = input('Take motors offline? y/n\n','s');
if strncmp(str,'y',1)
    delete(stage);
else
    return
end