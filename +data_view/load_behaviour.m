%% Manual define info
animal_id = 1;
rec_day = '2023-05-09';
rec_time = '1629';

% exp_info = general.AM_load_experiment(animal, rec_day, rec_time, bhv=true);

%% Save task behaviour 

animals = {'AP004','AP005'};

behaviour = struct;

% create matrix of times 
timestep = 0.01;
start_time = -2;
end_time = 2;
timevec = start_time:timestep:end_time;


behaviour.timestep = timestep;
behaviour.start_time = start_time;
behaviour.end_time = end_time;
behaviour.timevec = timevec;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    behaviour(animal_id).animal = animal;
    
    
    % find task days    
    protocol = 'stim_wheel_right*';
    experiments = general.find_experiments(animal,protocol);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;

        % save
        behaviour(animal_id).day{day_index} = day;
        experiment = experiments(day_index).experiment{end};
        
        % load experiment
        exp_info = general.AM_load_experiment(animal, day, experiment, bhv=true);
        
        % time vec around stim onsets
%         stimOn_times = exp_info.photodiode_times(2:end-1);
        stimOn_times = exp_info.photodiode_times;
        stimOn_times = stimOn_times(1:2:end);
        time_stimulus = stimOn_times+timevec;
        behaviour(animal_id).time_stimulus{day_index} = time_stimulus;
   
        % define t
        t =exp_info.timelite.timestamps;
        behaviour(animal_id).t{day_index} = t;

        % wheel position
        stim_wheel_position = interp1(t,exp_info.wheel_position,time_stimulus');
        behaviour(animal_id).stim_wheel_position{day_index} = stim_wheel_position;
        
        % stim aligned wheel move
        stim_wheel_move = interp1(t,+exp_info.wheel_move,time_stimulus');
        behaviour(animal_id).stim_wheel_move{day_index} = stim_wheel_move;
        
        % all moves
        tmp_move = [0; diff(exp_info.wheel_move)];
        all_move_on_frames = find(tmp_move==1);
        all_move_on_times = t(all_move_on_frames);
        behaviour(animal_id).all_move_on_times{day_index} = all_move_on_times;
                
        % all move offsets 
        all_move_off_frames = find(tmp_move==-1);
        all_move_off_times = t(all_move_off_frames);
        behaviour(animal_id).all_move_off_times{day_index} = all_move_off_times;
        
        % move after stim times
        stim_move_on_times = all_move_on_times(cell2mat(arrayfun(@(X) find(all_move_on_times>X,1,'first'),stimOn_times, 'UniformOutput', 0)));
        if length(stimOn_times)>length(stim_move_on_times)
            end_nans = nan(length(stimOn_times)-length(stim_move_on_times),1);
            stim_move_on_times = vertcat(stim_move_on_times, end_nans);
        end
        behaviour(animal_id).stim_move_on_times{day_index} = stim_move_on_times;
        
        
        % move offsets after stim times
        stim_move_off_times = all_move_off_times(cell2mat(arrayfun(@(X) find(all_move_on_times>X,1,'first'),stim_move_on_times, 'UniformOutput', 0)));
        behaviour(animal_id).stim_move_off_times{day_index} = stim_move_off_times;
        
        % reaction times
        reaction_times = stim_move_on_times - stimOn_times;
        behaviour(animal_id).reaction_times{day_index} = reaction_times;
        
    end
    
    disp(['Done with ' animal])
end

% end
disp('Done all')

% Save 
save('all_mice_behaviour.mat', 'behaviour', '-v7.3')
disp('Saved')

