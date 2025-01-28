save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andrada-Maria_Marica\save_stuff';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

workflow = {'stim_wheel_right*'};

stimwheel_pval = cell(length(animals), 12);
stimwheel_rxn_med = cell(length(animals), 12);
stimwheel_rxn_null_med = cell(length(animals), 12);
stimwheel_rnx_mean = cell(length(animals), 12);

bhv = struct;

for animal_idx=1:length(animals)
    animal = animals{animal_idx};
    recordings = plab.find_recordings(animal, [], workflow);

    for use_rec=1:length(recordings)

        rec_day = recordings(use_rec).day;
        rec_time = recordings(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        ap.load_recording

        [stimwheel_pval{animal_idx, use_rec},...
            stimwheel_rxn_med{animal_idx, use_rec}, ...
            stimwheel_rxn_null_med{animal_idx, use_rec}] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move);
        stimwheel_rxn_mean{animal_idx, use_rec} = mean(stim_to_move);

        learned_days{animal_idx, use_rec} = stimwheel_pval{animal_idx, use_rec} < 0.05;

        if use_rec==1
            learned_days{animal_idx, use_rec} = 0;
        end
    end

    learned_day = find([learned_days{animal_idx, :}], 1, 'first');
    recording_days = 1:length(recordings);
    days_from_learning{animal_idx} = recording_days - learned_day;

    bhv_days{animal_idx} = {recordings.day};

end

figure;
for animal_idx=1:length(animals)
    plot([stimwheel_pval{animal_idx,:}]);
    hold on
end
legend(animals)

bhv.animals = animals;
bhv.bhv_days = bhv_days;
bhv.stimwheel_pval = stimwheel_pval;
bhv.stimwheel_rxn_med = stimwheel_rxn_med;
bhv.stimwheel_rxn_null_med = stimwheel_rxn_null_med;
bhv.stimwheel_rxn_mean = stimwheel_rxn_mean;
bhv.learned_days = learned_days;
bhv.days_from_learning = days_from_learning;

save_name = fullfile(save_path, 'swr_bhv');
save(save_name, "bhv", "-v7.3");

clear all;