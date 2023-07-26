function out = AM_load_experiment(animal, rec_day, rec_time, load_parts)
% Loads data from experiments

% Define function arguments
arguments
    animal 
    rec_day 
    rec_time 
    load_parts.wf logical = false
    load_parts.bhv logical = false
    load_parts.ephys logical = false
    load_parts.mousecam logical = false
end

%% Load timelite and associated inputs

timelite_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'timelite.mat');

try
    timelite = load(timelite_fn);
    disp('Loading timelite...');
catch
    error(['No Timelite for: ' animal ' ' rec_day ' ' rec_time]);
end

% Set level for TTL threshold
ttl_thresh = 2;

% Flipper times
flipper_idx = strcmp({timelite.daq_info.channel_name}, 'flipper');
flipper_thresh = timelite.data(:,flipper_idx) >= ttl_thresh;
flipper_times = timelite.timestamps(find(diff(flipper_thresh) ~= 0) + 1);

% Mousecam times (note: all exposures, including those not recorded)
mousecam_idx = strcmp({timelite.daq_info.channel_name}, 'mouse_camera');
mousecam_thresh = timelite.data(:,mousecam_idx) >= ttl_thresh;
mousecam_expose_times = timelite.timestamps(find(diff(mousecam_thresh) == 1) + 1);

% Widefield times
widefield_idx = strcmp({timelite.daq_info.channel_name}, 'widefield_camera');
widefield_thresh = timelite.data(:,widefield_idx) >= ttl_thresh;
widefield_expose_times = timelite.timestamps(find(diff(widefield_thresh) == 1) + 1);

% Wheel position and velocity
timelite_wheel_idx = strcmp({timelite.daq_info.channel_name}, 'wheel');
wheel_position = timelite.data(:,timelite_wheel_idx);
[wheel_velocity,wheel_move] = plab.bhv.parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

% Screen on times
screen_idx = strcmp({timelite.daq_info.channel_name}, 'stim_screen');
screen_on = timelite.data(:,screen_idx) > ttl_thresh;

% Photodiode flips (interpolate from previous across screen flicker)
photodiode_thresh_level = 1; % low: has a relatively slow rise time
photodiode_idx = strcmp({timelite.daq_info.channel_name}, 'photodiode');
photodiode_thresh_screen_on = medfilt1(timelite.data(screen_on,photodiode_idx),3);
photodiode_thresh = interp1(timelite.timestamps(screen_on),photodiode_thresh_screen_on, ...
    timelite.timestamps,'previous','extrap') > photodiode_thresh_level;
photodiode_times = timelite.timestamps(find(diff(photodiode_thresh) ~= 0) + 1);
stimOn_times = photodiode_times(1:2:end);

% Reward times (if on past certain time, valve signal flips rapidly to
% avoid burnout - take the reward onset as flipping up and staying high for
% some length of samples)
reward_idx = strcmp({timelite.daq_info.channel_name}, 'reward_valve');
reward_thresh = timelite.data(:,reward_idx) >= ttl_thresh;
reward_on_pattern = [0,ones(1,3)]; % flip up, be consecutively high
reward_times = timelite.timestamps(strfind(reward_thresh',reward_on_pattern));

%% Load Bonsai
if load_parts.bhv == true
    try
        bonsai_file = 'bonsai_events.csv';
        bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai',bonsai_file);
        trial_events = AP_load_bonsai(bonsai_fn);
    catch me
        error(['No Bonsai for: ' animal ' ' rec_day ' ' rec_time]);
    end
end

% Sparse noise: get noise locations and times
try
    bonsai_file = 'NoiseLocations.bin';
    bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai',bonsai_file);
    fid = fopen(bonsai_fn);

    n_x_squares = trial_events.parameters.ScreenExtentX./trial_events.parameters.StimSize;
    n_y_squares = trial_events.parameters.ScreenExtentY./trial_events.parameters.StimSize;

    noise_locations = reshape(fread(fid),n_y_squares,n_x_squares,[]);
    fclose(fid);

    % Get stim times from photodiode (extrapolate: sparse noise photodiode
    % flips every N stim to give a more robust signal)
    photodiode_stim_idx = 1:trial_events.parameters.NthPhotodiodeFlip:size(noise_locations,3);
    % (check that the number of photodiode flips is expected)
    if length(photodiode_stim_idx) ~= length(photodiode_times)
        error('Sparse noise: mismatch photodiode times and stim number')
    end

    stim_times = interp1(photodiode_stim_idx,photodiode_times, ...
        1:size(noise_locations,3),'linear','extrap')';

catch me
end
%% Load mousecam
if load_parts.mousecam == true

mousecam_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

% Read mousecam header and get flipper times
mousecam_flipper_pin = 2;
mousecam_header = plab.mousecam.read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

% Check that timelite and mousecam have equal flipper flips
if length(flipper_times) ~= length(mousecam_flipper_times)
    warning('Flipper times not matched in timelite and mousecam')
    min_flipper_n = min(length(flipper_times),length(mousecam_flipper_times));
    %%% temporary?
    flipper_times = flipper_times(1:min_flipper_n);
    mousecam_flipper_times = mousecam_flipper_times(1:min_flipper_n);
end

% Get frame time after flips in timeline and mousecam
mousecam_postflips_idx_tl = arrayfun(@(x) ...
    find(mousecam_expose_times > flipper_times(x),1), ...
    1:length(flipper_times))';

mousecam_postflips_idx_cam = find(diff(mousecam_header.flipper) ~= 0) + 1;

% For sync: only use frames where flip happened in window before frame
% started (if flip happens during/close to exposure - camera pin state can
% be ambiguous)
mousecam_use_flips = ...
    (mousecam_expose_times(mousecam_postflips_idx_tl) - flipper_times) > 0.005 & ...
    (mousecam_expose_times(mousecam_postflips_idx_tl) - flipper_times) < 0.02;

use_flipframes = setdiff(1:length(flipper_times), ...
    find(diff(mousecam_postflips_idx_tl) ~= diff(mousecam_postflips_idx_cam)) + [0,1]);

% Get offset between frame index in timelite and mousecam
mousecam_idx_offset = unique( ...
    mousecam_postflips_idx_tl(mousecam_use_flips) - ...
    mousecam_postflips_idx_cam(mousecam_use_flips));

% If there's more than one offset value, something's misaligned
if length(mousecam_idx_offset) ~= 1
    error('Mousecam frames misaligned: >1 offset value')
end

% Get the corresponding timelite frame times for each mousecam frame
mousecam_tl_idx = (1:length(mousecam_header.timestamps)) + mousecam_idx_offset;
mousecam_times = mousecam_expose_times(mousecam_tl_idx);

end
%% Load widefield

if load_parts.wf == true

    % Load widefield data for all colors
    widefield_colors = {'blue','violet'};
    [wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
    for curr_wf = 1:length(widefield_colors)
        mean_image_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
            'widefield',sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
        svdU_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
            'widefield',sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
        svdV_fn = plab.locations.make_server_filename(animal,rec_day,rec_time, ...
            'widefield',sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

        wf_avg_all{curr_wf} = readNPY(mean_image_fn);
        wf_U_raw{curr_wf} = readNPY(svdU_fn);
        wf_V_raw{curr_wf} = readNPY(svdV_fn);
        % Assume colors go in order: dictated by Arduino
        wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
    end

    % Correct hemodynamics
    V_neuro_hemocorr = plab.wf.hemo_correct( ...
        wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
        wf_U_raw{2},wf_V_raw{2},wf_t_all{2});

    % Get DF/F
    wf_Vdf = plab.wf.svd_dff(wf_U_raw{1},V_neuro_hemocorr,wf_avg_all{1});

    % Deconvolve
    wf_framerate = mean(1./diff(wf_t_all{1}));
    wf_Vdf_deconv = AP_deconv_wf(wf_Vdf,[],wf_framerate);

    % Set final processed widefield variables
    wf_U = wf_U_raw{1};
    wf_V = wf_Vdf_deconv;
    wf_times = wf_t_all{1};
    wf_avg = wf_avg_all{1};

end

%% Load ephys

ephys_quality_control = false;

if load_parts.ephys == true

    ephys_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys','pykilosort');
    disp('Loading ephys...');

    % get probe folders
    probe_folders = dir(ephys_path);
    probe_folders = {probe_folders(3:end).name};
    probe_paths = strcat(ephys_path, filesep, probe_folders);

    for probe_idx=1:length(probe_paths)

        probe_path = probe_paths{probe_idx};

        % check contents are there and skip if not
        if length(dir(probe_path))<26
            continue
        end

        % These are the digital channels going into the FPGA
        flipper_sync_idx = 1;

        % Load phy sorting if it exists
        % (old = cluster_groups.csv, new = cluster_group.tsv because fuck me)
        cluster_filepattern = [probe_path filesep 'cluster_group*'];
        cluster_filedir = dir(cluster_filepattern);
        if ~isempty(cluster_filedir)
            cluster_filename = [probe_path filesep cluster_filedir.name];
            fid = fopen(cluster_filename);
            cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
            fclose(fid);
        end

        % Load sync/photodiode
        load(([probe_path filesep 'sync.mat']));

        % Read header information
        header_path = [probe_path filesep 'dat_params.txt'];
        header_fid = fopen(header_path);
        header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
        fclose(header_fid);

        header = struct;
        for i = 1:length(header_info{1})
            header.(header_info{1}{i}) = header_info{2}{i};
        end

        % Load spike data
        if isfield(header,'sample_rate')
            ephys_sample_rate = str2num(header.sample_rate);
        elseif isfield(header,'ap_sample_rate')
            ephys_sample_rate = str2num(header.ap_sample_rate);
        end
        spike_times = double(readNPY([probe_path filesep 'spike_times.npy']))./ephys_sample_rate;
        spike_templates_0idx = readNPY([probe_path filesep 'spike_templates.npy']);
        templates_whitened = readNPY([probe_path filesep 'templates.npy']);
        channel_positions = readNPY([probe_path filesep 'channel_positions.npy']);
        channel_map = readNPY([probe_path filesep 'channel_map.npy']);
        winv = readNPY([probe_path filesep 'whitening_mat_inv.npy']);
        template_amplitudes = readNPY([probe_path filesep 'amplitudes.npy']);

        % Default channel map/positions are from end: make from surface
        % (hardcode this: kilosort2 drops channels)
        max_depth = 3840;
        channel_positions(:,2) = max_depth - channel_positions(:,2);

        % Unwhiten templates
        templates = zeros(size(templates_whitened));
        for t = 1:size(templates_whitened,1)
            templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
        end

        % Get the waveform of all templates (channel with largest amplitude)
        [~,max_site] = max(max(abs(templates),[],2),[],3);
        templates_max = nan(size(templates,1),size(templates,2));
        for curr_template = 1:size(templates,1)
            templates_max(curr_template,:) = ...
                templates(curr_template,:,max_site(curr_template));
        end
        waveforms = templates_max;

        % Get depth of each template
        % (get min-max range for each channel)
        template_chan_amp = squeeze(range(templates,2));
        % (zero-out low amplitude channels)
        template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
        template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
        % (get center-of-mass on thresholded channel amplitudes)
        template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);

        % Get the depth of each spike (templates are zero-indexed)
        spike_depths = template_depths(spike_templates_0idx+1);

        % Get trough-to-peak time for each template
        templates_max_signfix = bsxfun(@times,templates_max, ...
            sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));

        [~,waveform_trough] = min(templates_max,[],2);
        [~,waveform_peak_rel] = arrayfun(@(x) ...
            max(templates_max(x,waveform_trough(x):end),[],2), ...
            transpose(1:size(templates_max,1)));
        waveform_peak = waveform_peak_rel + waveform_trough;

        templateDuration = waveform_peak - waveform_trough;
        templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;

        % Get sync points for alignment

        % Get index of this experiment within day
        % (to determine which flipper boundaries match this recording)
        experiment_idx = find(strcmp(ap.find_recordings(animal,rec_day).protocol,rec_time));

        if exist('flipper_times','var')
            % (if flipper, use that)
            % (at least one experiment the acqLive connection to ephys was bad
            % so it was delayed - ideally check consistency since it's
            % redundant)
            bad_flipper = false;

            % Get flipper experiment differences by long delays
            % (note: this is absolute difference, if recording stopped and
            % started then the clock starts over again, although I thought it
            % wasn't supposed to when I grab the concatenated sync, so
            % something might be wrong)
            flip_diff_thresh = 10; % time between flips to define experiment gap (s)
            flipper_expt_idx = [1;find(abs(diff(sync(flipper_sync_idx).timestamps)) > ...
                flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];

            flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
                flipper_expt_idx(experiment_idx):flipper_expt_idx(experiment_idx+1)-1);

            % Pick flipper times to use for alignment
            if length(flipper_flip_times_ephys) == length(flipper_times)
                % If same number of flips in ephys/timeline, use all
                sync_timeline = flipper_times;
                sync_ephys = flipper_flip_times_ephys;
            elseif length(flipper_flip_times_ephys) ~= length(flipper_times)
                % If different number of flips in ephys/timeline, best
                % contiguous set via xcorr of diff
                warning([animal ' ' rec_day ':Flipper flip times different in timeline/ephys'])
                warning(['The fix for this is probably not robust: always check'])
                [flipper_xcorr,flipper_lags] = ...
                    xcorr(diff(flipper_times),diff(flipper_flip_times_ephys));
                [~,flipper_lag_idx] = max(flipper_xcorr);
                flipper_lag = flipper_lags(flipper_lag_idx);
                % (at the moment, assuming only dropped from ephys)
                sync_ephys = flipper_flip_times_ephys;
                sync_timeline = flipper_times(flipper_lag+1: ...
                    flipper_lag+1:flipper_lag+length(flipper_flip_times_ephys));
            end

        else
            bad_flipper = true;
        end

        if bad_flipper
            % (if no flipper or flipper problem, use acqLive)

            % Get acqLive times for current experiment
            experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
            experiment_ephys_stops = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 0);
            acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_idx), ...
                experiment_ephys_stops(experiment_idx)];

            sync_timeline = acqLive_timeline;
            sync_ephys = acqlive_ephys_currexpt;

            % Check that the experiment time is the same within threshold
            % (it should be almost exactly the same)
            if abs(diff(acqLive_timeline) - diff(acqlive_ephys_currexpt)) > 1
                error([animal ' ' day ': acqLive duration different in timeline and ephys']);
            end
        end

        % Get spike times in timeline time
        spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');

        disp(['Loaded probe ' num2str(probe_idx)]);
    end
    
    %%%%%%%%% BOMBCELL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% WORKING ON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qMetricsExist = 0;
    qMetrics_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys','qMetrics');
    if isfolder(qMetrics_path)
        qMetricsExist = 0;
        qMetrics_probe_folders = dir(qMetrics_path);
        qMetrics_probe_folders = {qMetrics_probe_folders(3:end).name};
        qMetrics_probe_paths = strcat(qMetrics_path, filesep, qMetrics_probe_folders);
        
        for probe_idx=1:length(qMetrics_probe_paths)

            qMetrics_probe_path = qMetrics_probe_paths{probe_idx};

            % Load Julie's quality metrics
            qMetricsExist = ~isempty(dir(fullfile(qMetrics_probe_path, 'qMetric*.mat'))) || ~isempty(dir(fullfile(qMetrics_probe_path, 'templates._bc_qMetrics.parquet')));
            if qMetricsExist
                [param, qMetric] = bc_loadSavedMetrics(qMetrics_probe_path);
                unitType = bc_getQualityUnitType(param, qMetric);
                disp(['Loaded qMetrics for probe ' num2str(probe_idx)]);

%                 % Define good units from labels
%                 good_templates_idx = find(unitType == 1 | unitType == 2);
%                 good_templates = ismember(0:size(templates,1)-1,good_templates_idx);       
            end
        end
    end


    %%%%% PHY SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get "good" templates from labels if quality control selected and
    % manual labels exist
    if ephys_quality_control && exist('cluster_groups','var') % && qMetricsExist == 0
        % If there's a manual classification
        if verbose; disp('Keeping manually labelled good units...'); end
        
        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled']);
        end
        
%         % Define good units from labels
%         good_templates_idx = uint32(cluster_groups{1}( ...
%             strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
%         good_templates = ismember(0:size(templates,1)-1,good_templates_idx);        
    else
        % If no cluster groups at all, keep all
        warning([animal ' ' rec_day ' - no ephys quality control']);
        disp('No ephys quality control, keeping all and re-indexing'); 
%         good_templates_idx = unique(spike_templates_0idx);
%         good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    end
    
%     % Throw out all non-good template data
%     templates = templates(good_templates,:,:);
%     template_depths = template_depths(good_templates);
%     waveforms = waveforms(good_templates,:);
%     templateDuration = templateDuration(good_templates);
%     templateDuration_us = templateDuration_us(good_templates);
%     
%     % Throw out all non-good spike data
%     good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
%     spike_times = spike_times(good_spike_idx);
%     spike_templates_0idx = spike_templates_0idx(good_spike_idx);
%     template_amplitudes = template_amplitudes(good_spike_idx);
%     spike_depths = spike_depths(good_spike_idx);
%     spike_times_timeline = spike_times_timeline(good_spike_idx);
%     
%     % Rename the spike templates according to the remaining templates
%     % (and make 1-indexed from 0-indexed)
%     new_spike_idx = nan(max(spike_templates_0idx)+1,1);
%     new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
%     spike_templates = new_spike_idx(spike_templates_0idx+1);
       
end
%% Experiment scroller (check alignment)

% AP_expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2},mousecam_fn,mousecam_times)


% AP_expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)

%%  Return everything in the workspace
workspace_vars = who;
out = struct;
for n_var = 1:length(workspace_vars)
    out.(workspace_vars{n_var}) = eval(workspace_vars{n_var});
end

end
