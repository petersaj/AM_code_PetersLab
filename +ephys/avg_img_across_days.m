
workflow = {'lcr_passive'};

% 'AM014', 'AM015',

animals = {'AP023'} %, 'AM027'};

% save_fig_path = 'C:\Users\amarica\Documents\Lab stuff\Random figs\Long_str_stuff\Passive\Probe_position';

for animal_idx=1:length(animals)
    animal = animals{animal_idx};

    recordings = plab.find_recordings(animal, [], workflow);

    fig_mean_img = figure('Position', get(0, 'Screensize'));
    mean_img = tiledlayout('flow');
    for use_rec = 4:length(recordings)
%         use_rec = length(recordings)-1;
        
        rec_day = recordings(use_rec).day;
        rec_time = recordings(use_rec).recording{end};
        
        verbose = true;
        load_parts.widefield = true;
        
        ap.load_recording
    
        violet_idx = strfind(mean_image_fn, 'violet');
        mean_image_fn_blue = [mean_image_fn(1:violet_idx-1) 'blue.npy'];
        mean_image_blue = readNPY(mean_image_fn_blue);

        nexttile
        imagesc(mean_image_blue)
        axis image
        colormap('gray')

        %     ap.wf_draw('grid', 'y', true)
        %     ap.wf_draw('point', [0.6, 0.7], true)

        clim([500, 20000])
    end
    mean_img.Title.String = [animal ' avg img across days'];
    mean_img.Title.FontSize = 18;
    mean_img.Title.FontWeight = "bold";

    % save fig
%     fig_mean_img_name = [animal '_18_19_avg_blue_img_across_days.tif'];
%     fig_mean_img_path = fullfile(save_fig_path, fig_mean_img_name);
%     saveas(fig_mean_img, fig_mean_img_path);
end