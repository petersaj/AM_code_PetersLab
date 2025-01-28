
animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

for animal_idx=1:length(animals)
    animal = animals{animal_idx};
    disp(['Start ' animal])
    recordings = plab.find_recordings(animal, []);

    ephys_recordings = find([recordings.ephys] == 1);

    for rec_idx=ephys_recordings

        day = recordings(rec_idx).day;
        %% Get paths and filenames

        ephys_path = plab.locations.filename('server',animal,day,[],'ephys');

        ephys_exists = exist(ephys_path,'dir');

        if ~ephys_exists
            error([animal ' ' day ': No ephys data found']);
        end

        save_paths = {[ephys_path filesep 'kilosort4']};
        data_paths = {ephys_path};

        % Get ephys recording paths
        % probe_n = multiple probes simultaneously
        % site_n = multiple sites recorded in serial
        data_path_dir = ...
            [dir(fullfile(data_paths{1}, 'probe_*')), ...
            dir(fullfile(data_paths{1}, 'site_*'))];
        if ~isempty(data_path_dir)
            data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir.name},'uni',false);
            save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir.name},'uni',false);
        end

        %% Run kilosort on all datasets

        for curr_data = 1:length(data_paths)

            % Get experiments
            % (multiple if preview off/on, default process separately)
            curr_data_path = data_paths{curr_data};
            ephys_exp_paths = dir([curr_data_path filesep 'experiment*']);

            for curr_exp = 1:length(ephys_exp_paths)

                curr_exp_path = fullfile(ephys_exp_paths(curr_exp).folder, ...
                    ephys_exp_paths(curr_exp).name);

                % Update save path with experiment (only if more than one, legacy)
                if length(ephys_exp_paths) == 1
                    curr_save_path = save_paths{curr_data};
                elseif length(ephys_exp_paths) > 1
                    curr_save_path = fullfile(save_paths{curr_data},ephys_exp_paths(curr_exp).name);
                end

                % Find Open Ephys recordings
                % (multiple if record off/on, default concatenate processing)
                experiment_dir = dir(fullfile(curr_exp_path,'recording*'));

                % Get Open Ephys filename(s)
                ap_data_dir = cellfun(@(data_path,data_fn) ...
                    dir(fullfile(data_path,data_fn,'continuous','*-AP','continuous.dat')), ...
                    {experiment_dir.folder},{experiment_dir.name},'uni',false);
                ap_data_filenames = cellfun(@(data_dir) ...
                    fullfile(data_dir.folder,data_dir.name),ap_data_dir,'uni',false);

                %% Get metadata filename (for sample rate: just use first file)
                ephys_meta_dir = dir(fullfile(experiment_dir(1).folder,experiment_dir(1).name,'**','*.oebin'));
                ephys_meta_fn = fullfile(ephys_meta_dir.folder,ephys_meta_dir.name);

                %% copy old folder
                qMetrics_path = fullfile(curr_save_path, 'qMetrics');
                old_qMetrics_path = fullfile(curr_save_path, 'old_qMetrics');
                mkdir(old_qMetrics_path)
                try
                    copyfile(qMetrics_path, old_qMetrics_path);
                    disp('Folder copied successfully!');
                catch ME
                    fprintf('Error copying folder: %s\n', ME.message);
                    fprintf('\n For: %s %s\n',animal,day);
                end
                %% Run bombcell

                % Run bombcell
                kilosort_version = 4;
                %%% MAKE SURE this has rerun = 1 %%%%%%%%%%%%%%%%%%%%%%%
                ap.run_bombcell(ap_data_filenames{1},curr_save_path,ephys_meta_fn,kilosort_version);

            end
        end

        %% Print end message
        fprintf('\nDone bombcell: %s %s\n',animal,day);
    end

end

