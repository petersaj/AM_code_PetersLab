%% ~~~~~~~~ ALIGN WIDEFIELD

%% Create animal alignment (animal average VFS to master VFS)
% NOTE: need day alignment and retinotopy

% already aligned, 14 and 15 don't have retinotopy
% animals = ["AM011"; "AM012"; "AM014"; "AM015"; "AM016"; ...
%     "AM017"; "AM018"; "AM019"; "AM021"; "AM022"; "AM026"];
animals = {'AP023', 'AP025'};
for animal_idx=1:length(animals)
    animal = animals{animal_idx};

    % Create across-day alignments
    plab.wf.wf_align([],animal,[],'new_days');

%     % Get and save VFS maps for animal
    plab.wf.retinotopy_vfs_batch(animal);

    % Create across-animal alignments
    plab.wf.wf_align([],animal,[],'new_animal');
end


%% View aligned days

animal = 'AM014';

recordings = plab.find_recordings(animal);
wf_days_idx = cellfun(@(x) any(x),{recordings.widefield});
wf_recordings = recordings(wf_days_idx);

avg_im_aligned = cell(size(wf_recordings));
for curr_day = 1:length(wf_recordings)
    day = wf_recordings(curr_day).day;

    img_path = plab.locations.filename('server', ...
        animal,day,[],'widefield');

    avg_im_n = readNPY([img_path filesep 'meanImage_blue.npy']);
    avg_im_h = readNPY([img_path filesep 'meanImage_violet.npy']);

%         % (to concatenate)
%         avg_im_aligned{curr_day} = [AM_wf_align(avg_im_n,animal,day), ...
%             AM_wf_align(avg_im_h,animal,day)];

    % (blue only)
    avg_im_aligned{curr_day} = plab.wf.wf_align(avg_im_n,animal,day);
end

% Plot average
c = prctile(reshape([avg_im_aligned{:}],[],1),[0,99.9]);
AP_imscroll(cat(3,avg_im_aligned{:}),{wf_recordings.day});
caxis(c);
axis image;
set(gcf,'Name',animal);


%% TO DRAW STUFF

animal = 'AM021'

recordings = plab.find_recordings(animal);
wf_days_idx = cellfun(@(x) any(x),{recordings.widefield});
wf_recordings = recordings(wf_days_idx);

curr_day = length(wf_recordings);
day = wf_recordings(curr_day).day;

img_path = plab.locations.filename('server', ...
    animal,day,[],'widefield');

avg_im_n = readNPY([img_path filesep 'meanImage_blue.npy']);
avg_im_h = readNPY([img_path filesep 'meanImage_violet.npy']);

% (blue only)
avg_im_aligned = AM_wf_align(avg_im_n,animal,day);

figure
imagesc(avg_im_aligned)
axis image
axis off
colormap('gray')

% ap.wf_draw('ccf','r')
ap.wf_draw('point', [0.5, 0.5])
ap.wf_draw('grid', 'y')

clim([5000, 25000])



