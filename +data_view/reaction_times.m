%% Load and plot reaction times

%% all mice
load('all_mice_behaviour.mat');

animals = {'AP004','AP005'};

% Histograms per day
for animal_id=1:length(animals)
    figure;
    tiledlayout("flow");
    animal = animals{animal_id};
    all_training_days = behaviour(animal_id).day;

    % get reaction times 
    for day_idx = 1:length(all_training_days)
        day = all_training_days{day_idx};
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        % plot
        nexttile;
        histogram(reaction_times,'BinWidth', 0.05, 'BinLimits', [0 1], 'Normalization', 'probability')
        hold on;
        xlabel('Reaction time')
        ylabel('Probability')
        title([animal ' ' day])
    end
end

% Summary RT plot
figure;
for animal_id=1:length(animals)
    all_training_days = behaviour(animal_id).day;

    % get reaction times percentage that falls between 100-200 ms
    reaction_times_percentage = nan(1,length(all_training_days));
    for day_idx = 1:length(all_training_days)
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        reaction_times_percentage(day_idx) = sum(discretize(reaction_times,[0.1 0.2])==1)/length(reaction_times);
    end
    % make plot
    plot(reaction_times_percentage,'o-')
    hold on;
end

%title('Reaction times for original task')
xlabel('Training day')
ylabel('% of fast reaction times')
legend(animals, 'Location','northwest')
AM_prettyfig



