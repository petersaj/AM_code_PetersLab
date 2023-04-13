% Loads data from experiments
%
% Settings:
% (what to load)
% load_parts.(cam/imaging/ephys)
%
% (ephys)
% site (if multiple probes)
% recording (if multiple on/off recordings within probe)

% Turn warnings off
warning on

%% Display progress or not
if ~exist('verbose','var')
    verbose = false;
end

%% Define what to load

% Site (multiple probes) is optional
if ~exist('site','var')
    site = [];
end

% Recording (multiple recordings on one probe) is optional
if ~exist('recording','var')
    recording = [];
end

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts','var')
    load_parts.cam = true;
    load_parts.imaging = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'cam')
        load_parts.cam = false;
    end
    if ~isfield(load_parts,'imaging')
        load_parts.imaging = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end


%% Manual define info

animal = 'new_imagetest';
day = '2023-04-13';
experiment = '1815';

verbose = true;

%% Load timeline and associated inputs

timelite_filename = plab.locations.make_server_filename(animal,day,experiment,'timelite.mat');

try
    load(timelite_filename);
    if verbose
        disp('Loading timelite...'); 
    end
catch
     error(['No Timelite for: ' animal ' ' day ' ' experiment]);
end

% Get widefield camera times
cam_name = 'widefield_camera';
timeline_cam_idx = strcmp({daq_info.channel_name}, cam_name);
wf_strobe = data(:,timeline_cam_idx);
cam_expose_times = timestamps(find(diff(wf_strobe>2)==1));

% Get flipper
timelite_flipper = data(:,strcmp({daq_info.channel_name}, 'flipper'));
flipper_trace = timelite_flipper>2;

% mousecam stuff
mousecam_header_fn = "P:\Data\new_imagetest\2023-04-13\Protocol_1815\mousecam_header.bin";
flipper_pin = 2;
mousecam_header = read_mousecam_header(mousecam_header_fn, flipper_pin);

good_mousecam_frames = find(mousecam_header.flipper>0, 1):length(mousecam_header.flipper);

%%%%%%%%%%%%%% get times ?????? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
facecamera_times = interp1(double(flipper_trace), ...
    timestamps, ...
    double(mousecam_header.flipper(good_mousecam_frames)));

% plot 
figure;
plot(mousecam_header.timestamps-mousecam_header.timestamps(1), mousecam_header.flipper);
hold on;
plot(timestamps, flipper_trace+1);


%%%%%%%%%%%%%%%%%%% NOT INCLUDED YET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Get acqLive signal
%     acqLive_name = 'acqLive';
%     acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
%     thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
%     acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
%     acqLive_timeline = Timeline.rawDAQTimestamps( ...
%         [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);
    
    % Get wheel position and velocity
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    % (this is a very strange hack to overcome a problem in the rotary
    % encoder that's known in the lab and was put on the wiki)
    wheel_position = Timeline.rawDAQData(:,rotaryEncoder_idx);
    wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;    
    [wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);
    
    % Get whether stim was flickering
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
    end
    
    % Get photodiode flips (compensate for screen flicker)
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    % (define stim screen on from photodiode - sometimes sample-length
    % offset maybe because of backlight onset delay)
    stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.2;
    stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
    photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
    % (medfilt because photodiode can be intermediate value when backlight
    % coming on)
    % (OLD: this worked fine until photodiode bug: sometimes gray)
%     photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
%         photodiode_idx),3) > photodiode_thresh;
%     photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
%         (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
%     photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
    % (NEW: accomodating photodiode bug flipping sometimes to gray)
    photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx),3);
    photodiode_diff_thresh = range(Timeline.rawDAQData(:,photodiode_idx))*0.2;
    photodiode_diff_t = 50; % time (in ms) to get delayed differential
    photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
    photodiode_diff_filt = [1,zeros(1,photodiode_diff_samples),-1];
    photodiode_trace_diff = abs(conv(photodiode_trace_medfilt,photodiode_diff_filt,'valid')) > ...
        photodiode_diff_thresh;
    photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
        photodiode_trace_diff(2:end))+ photodiode_diff_samples + 1;
    photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

    % Get flipper signal (this was added late, might not be present)
    flipper_name = 'flipper';
    flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
    flipper_thresh = 2; % TTL threshold
    flipper_trace = Timeline.rawDAQData(:,flipper_idx) > flipper_thresh;
    flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
        (flipper_trace(1:end-1) & ~flipper_trace(2:end)))+1;
    flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';
    
end