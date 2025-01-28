% workflow = 'lcr_passive';
% animals = {'AP008'};
% rerun = 1;
% 
% for animal_idx=1:length(animals)
%     animal = animals{animal_idx};
% 
%     experiments = ap.find_recordings(animal, [], workflow)
%     days = {experiments.day};
%     for day_idx=1:length(days)
%         rec_day =  days{day_idx};
%         ephys.run_bombcell(animal,rec_day, rerun)
%     end
% end
% 


workflow = 'lcr_passive';
animals = {'AP009'};
% rerun = 1;

for animal_idx=1:length(animals)
    animal = animals{animal_idx};

    experiments = ap.find_recordings(animal, [], workflow)
    days = {experiments.day};
    for day_idx=4:length(days)
        rec_day =  days{day_idx};
        ap.run_bombcell(animal,rec_day)
    end
end


