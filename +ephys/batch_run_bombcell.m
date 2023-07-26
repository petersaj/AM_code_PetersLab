protocol = 'lcr_passive';
animals = {'AP009'};
rerun = 1;

for animal_idx=1:length(animals)
    animal = animals{animal_idx};

    experiments = general.find_experiments(animal, protocol)
    days = {experiments.day};
    for day_idx=1:length(days)
        rec_day =  days{day_idx};
        ephys.run_bombcell(animal,rec_day, rerun)
    end
end

