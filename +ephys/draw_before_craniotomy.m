animal = 'AP025'

workflow = {'lcr_passive'};


recordings = plab.find_recordings(animal, [], workflow);

use_rec = length(recordings);

rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).recording{end};

verbose = true;

ap.load_recording

% a = readNPY("P:\Data\AM028\2024-08-12\widefield\meanImage_blue.npy");
figure
imagesc(wf_avg)
axis image
colormap('gray')

ap.wf_draw('grid', 'y')
ap.wf_draw('point', [0.6, 0.7])

clim([0, 20000])
axis off