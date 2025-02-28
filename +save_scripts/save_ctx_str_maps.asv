%% save path and animals in dataset

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\long_str_ctx_data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

for animal_idx=1:length(animals)
    animal = animals{animal_idx};
    disp(['Start ' animal])
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);
    bhv_days = {train_rec_passive.day};
%     wf_days =  bhv_days([train_rec_passive.widefield]);
%     ephys_days =  bhv_days([train_rec_passive.ephys]);
%     
    ctx_maps_to_str = table;
    for use_rec=1:length(train_rec_passive)
        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        ap.load_recording;

%         this_day_from_learning = days_from_learning(ismember(bhv_days, rec_day));

  
        %% Get striatum MUA in bins of set lengths
        % Set MUA bin length (microns)
        mua_length = 200;

        AP_longstriatum_find_striatum_depth
        mua_length = 200;
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(end);
        try
            depth_group = discretize(spike_depths,depth_group_edges);
        catch ME
            warning(['Undefined str depth ' animal ' ' rec_day])
            %             depth_group = nan;
            ctx_maps_to_str.animal(use_rec) = {animal};
            ctx_maps_to_str.rec_day(use_rec) = {rec_day};
            continue
        end

        % Get time bins corresponding to widefield frame exposures
        % (skip the beginning and end of the recording to avoid artifacts)
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        % Bin spikes in depth and time
        binned_spikes = zeros(max(depth_group),length(time_bins)-1);
        for curr_depth = 1:max(depth_group)
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        %% Regress ctx fluorescence to striatal MUA

        % Normalize MUA by standard devation
        % (to have kernels in comparable units)
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

        % Set parameters for regression
        % (widefield components to use - too many can overload memory or add noise,
        % too few gives a less detailed map)
        use_svs = 1:500;
        % (frame lags to use: for maps rather than prediction, can just do zero-lag
        % regression with no temporal component)
        kernel_t = [0,0];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        % (ridge regression normalizing factor lambda: too low gives just noise,
        % higher removes map detail)
        lambda = 5;
        % (whether to z-score regressors or signals: not here, since everything
        % already normalized)
        zs = [false,false];
        % (cross validation: not needed for maps, only for prediction)
        cvfold = 1;
        % (whether to use (yes) and return as ouput (no) a constant in regression)
        use_constant = true;
        return_constant = false;

        % Interpolate the widefield to the center of each time bin to match the
        % MUA timepoints
        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Do regression from cortex to spikes
        [cortex_kernel,predicted_spikes,explained_var] = ...
            ap.regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % (note: cortex_kernel is components x time x MUA - so in this example is
        % 500 used components x 1 time point x 12 MUA signals)

%         % Convert ctx kernal V into pixels, plot
         cortex_kernel_px = squeeze(plab.wf.svd2px(wf_U(:,:,use_svs),cortex_kernel));
%         ap.imscroll(cortex_kernel_px);
%         axis image;
%         clim(max(abs(clim)).*[-1,1]*0.7);
%         ap.wf_draw('ccf','k');
%         colormap(ap.colormap('PWG'));


        %% save everything
        ctx_maps_to_str.animal(use_rec) = {animal};
        ctx_maps_to_str.rec_day(use_rec) = {rec_day};

        ctx_maps_to_str.depth_group_edges(use_rec) = {depth_group_edges};
        ctx_maps_to_str.cortex_kernel_px(use_rec) = {cortex_kernel_px}; 
        ctx_maps_to_str.explained_var(use_rec) = {explained_var.total};

        disp(['Done day ' num2str(use_rec)])
        
    end
    all_ctx_maps_to_str_cell{animal_idx} = ctx_maps_to_str;
    disp(['Done ' animal])
end

all_ctx_maps_to_str = vertcat(all_ctx_maps_to_str_cell{:});

save_name = fullfile(save_path, 'ctx_maps_to_str');
save(save_name, "all_ctx_maps_to_str", "-v7.3");

