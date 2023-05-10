function experiments = find_experiments(animal, protocol)
% Finds experiments for specific animal and protocol

% find all days for this animal
animal_data_path = fullfile(plab.locations.server_data_path, animal);
animal_data_dir = dir(animal_data_path);
day_idx = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{animal_data_dir.name}) &...
    [animal_data_dir.isdir];
days = {animal_data_dir(day_idx).name};
days_pathnames = cellfun(@(x) [animal_data_path filesep x],days,'uni',false);

protocol_expts = cell(size(days));
% imaging_expts = cell(size(days));
% ephys_expts = cell(size(days));

for day_idx = 1:length(days)

    day = days{day_idx};

    % Find all protocols of that day
    curr_day_dir = dir(days_pathnames{day_idx});

    protocol_format = ['Protocol_(?<time>\d+)'];
    protocol_temp =  regexp({curr_day_dir.name}, protocol_format, 'names');
    protocol_temp = [protocol_temp{1,~cellfun(@isempty, protocol_temp)}];
    protocol_time = {protocol_temp.time};

    protocol_nums = length(protocol_time);

    % If looking for specific protocol, find amongst days's experiments
    if ~isempty(protocol)
        use_exp = false(1, protocol_nums);
        for curr_exp = 1:protocol_nums
            % Check for bonsai
            protocol_folder = plab.locations.make_server_filename(animal,day,protocol_time{curr_exp},'bonsai');
            is_good_protocol = ~isempty(dir(fullfile(protocol_folder,protocol)));
            if is_good_protocol
                use_exp(curr_exp) = 1;
            else
                continue
            end
        end
    else
        use_exp = true(1, protocol_nums);
    end

    protocol_expts{day_idx} = protocol_time(use_exp);
end


use_days = ~cellfun(@isempty,protocol_expts);

experiments = struct('day',cell(sum(use_days),1),'experiment',cell(sum(use_days),1));

[experiments.day] = deal(days{use_days});
[experiments.experiment] = deal(protocol_expts{use_days});
end