%% take out noisy units

% across all non-zero channels 
new_raw_waveforms = permute(rawWaveforms.average, [1 3 2]);
all_unit_non_zero_corr = nan(size(new_raw_waveforms,1), 1);
for unit_idx=1:size(new_raw_waveforms, 1)

    % get non-zero channels in template waveforms
    temp_template_waveforms = squeeze(templateWaveforms(unit_idx,:,:));

    [template_max_amplitude, template_max_channel] = max(max(abs(templateWaveforms(unit_idx, :, :)), [], 2), [], 3);

    zero_channel_idx = arrayfun(@(X) all(temp_template_waveforms(:,X) < 0.05* template_max_amplitude...
        & temp_template_waveforms(:,X) > -0.05* template_max_amplitude), [1:param.nChannels]);

    non_zero_template_waveforms = temp_template_waveforms(:, ~zero_channel_idx);

    % concatenate for comparison
    vector_non_zero_template_waveforms = non_zero_template_waveforms(:);

    % get raw waveforms on the same channels
    temp_raw_waveforms = squeeze(new_raw_waveforms(unit_idx,:,:));
    non_zero_raw_waveforms = temp_raw_waveforms(:, ~zero_channel_idx);
    vector_non_zero_raw_waveforms = non_zero_raw_waveforms(:);

    all_unit_non_zero_corr(unit_idx) = corr(vector_non_zero_raw_waveforms, vector_non_zero_template_waveforms);
end

% figure;
% histogram(all_unit_non_zero_corr,50)
% 
% figure;
% imagesc(squeeze(templateWaveforms(unit_idx, :, :)))

cutoff_all_unit_non_zero_corr = 0.6;
potential_noise_non_zero = find(all_unit_non_zero_corr<cutoff_all_unit_non_zero_corr)';


% compare max channels on raw and template
all_units_corr_max_chan = nan(size(new_raw_waveforms, 1), 1);

for unit_idx=1:size(templateWaveforms,1)
    % find channel with max amplitude
    [template_max_amplitude, template_max_channel] = max(max(abs(templateWaveforms(unit_idx, :, :)), [], 2), [], 3);
    [raw_max_amplitude, raw_max_channel] = max(max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2), [], 3);

    % get corr with max raw channel 
    raw_wvf = squeeze(new_raw_waveforms(unit_idx, :, raw_max_channel));
    template_wvf = squeeze(templateWaveforms(unit_idx, :, template_max_channel));
    all_units_corr_max_chan(unit_idx) = corr(raw_wvf', template_wvf');
end
% 
% figure;
% histogram(all_units_corr_max_chan,50)

cutoff_all_units_corr_max_chan = 0;
potential_noise_max_chan = find(all_units_corr_max_chan<cutoff_all_units_corr_max_chan)';

% compare max channels on raw and template
all_units_amp_decay = nan(size(new_raw_waveforms, 1), 1);

for unit_idx=1:size(new_raw_waveforms,1)
    % find channel with max amplitude
    [raw_max_amplitude, raw_max_channel] = max(max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2), [], 3);

    raw_max_ampl_all_channels = max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2);
        
    all_units_amp_decay(unit_idx) = sum(raw_max_ampl_all_channels>0.7*raw_max_amplitude);
end

cutoff_all_units_amp_decay = 10;
potential_noise_amp_decay = find(all_units_amp_decay>cutoff_all_units_amp_decay)';

% figure; histogram(all_units_amp_decay, 100)

% get rid off all noisy units
potential_noise = unique([potential_noise_max_chan potential_noise_non_zero potential_noise_amp_decay]);

am_bc_units = unitType;
am_bc_units(potential_noise) = 0;

extra_qmetric = struct;
extra_qmetric.all_unit_non_zero_corr = all_unit_non_zero_corr;
extra_qmetric.cutoff_all_unit_non_zero_corr = cutoff_all_unit_non_zero_corr;
extra_qmetric.all_units_corr_max_chan = all_units_corr_max_chan;
extra_qmetric.cutoff_all_units_corr_max_chan = cutoff_all_units_corr_max_chan;
extra_qmetric.all_units_amp_decay = all_units_amp_decay;
extra_qmetric.cutoff_all_units_amp_decay = cutoff_all_units_amp_decay';

% plots for extra metrics
figure;
potential_noise_plot = tiledlayout('flow');
title(potential_noise_plot, 'Potential noise')
nexttile
histogram(extra_qmetric.all_units_corr_max_chan,50)
hold on;
xline(extra_qmetric.cutoff_all_units_corr_max_chan, 'r');
title('Corr max channels')
nexttile
histogram(extra_qmetric.all_unit_non_zero_corr,50)
hold on;
xline(extra_qmetric.cutoff_all_unit_non_zero_corr, 'r');
title('Corr non-zero channels')
nexttile
histogram(extra_qmetric.all_units_amp_decay,50)
hold on;
xline(extra_qmetric.cutoff_all_units_amp_decay, 'r');
title('Amp decay')
