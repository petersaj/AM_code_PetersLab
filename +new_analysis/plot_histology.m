%% Grab and plot histology pictures (SMZ)
unique_animals = unique(bhv.animal);

for animal_idx=1:length(unique_animals)

    animal = unique_animals{animal_idx};
    % Just load all images
    histology_path = plab.locations.filename('server',animal,[],[],'histology');
    histology_dir = dir(fullfile(histology_path, 'processed_smz', '*.tif'));
    if isempty(histology_dir)
        histology_dir = dir(fullfile(histology_path, '*.tif'));
    end
    histology_filenames = cellfun(@(path,name) fullfile(path,name), ...
        {histology_dir.folder},{histology_dir.name},'uni',false);
    [~,sort_idx] = natsortfiles(histology_filenames);
    histology_im = cell(size(histology_dir));
    for curr_im = 1:length(sort_idx)
        histology_im{curr_im} = tiffreadVolume( ...
            histology_filenames{sort_idx(curr_im)});
    end
    n_chan = size(histology_im{1},3);
    % % Plot channel montage separately
    % figure('Name',animal); tiledlayout(1,n_chan);
    % for curr_chan = 1:n_chan
    %     nexttile;
    %     m = montage(cellfun(@(im) im(:,:,curr_chan),histology_im,'uni',false));
    %     clim(mean(m.CData,[1,2])*[0.25,4]);
    % end
    % Grab image montage and display as RGB
    im_montage = uint16([]);
    chan_cols = [0,1,0;1,0,0];
    for curr_chan = 1:n_chan
        h = figure;
        m = montage(cellfun(@(im) im(:,:,curr_chan),histology_im,'uni',false));
        im_montage = cat(3,im_montage,m.CData);
        close(h);
    end
    m_clim = arrayfun(@(chan) double(prctile(im_montage(:,:,chan),[20,90],'all')).*[1;1.5],1:n_chan,'uni',false);
    im_montage_rgb = min(sum(cell2mat(arrayfun(@(chan) ...
        mat2gray(im_montage(:,:,chan),double(m_clim{chan})).* ...
        permute(chan_cols(chan,:),[1,3,2]), ...
        permute(1:n_chan,[1,3,4,2]),'uni',false)),4),1);
    figure;image(im_montage_rgb);axis image off;
    title(['Animal ' num2str(animal_idx)]);
end