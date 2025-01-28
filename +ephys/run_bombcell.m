function run_bombcell(animal,rec_day, rerun)
%% run bombcell per recording
%% set paths
% find ephys recording path
% ephys_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys');
ephys_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys');
ephys_rec = dir(fullfile(ephys_path, 'experiment*', 'recording*'));

if ~isempty(ephys_rec)
    % get meta/oebin file
    ephys_meta_dir = dir(fullfile(ephys_rec.folder, ephys_rec.name, '*.oebin'));

    % get raw data
    ephys_raw_dir = dir(fullfile(ephys_rec.folder, ephys_rec.name, 'continuous', '*AP', '*.dat'));

    % kilosort folder
    ephys_ks_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys','pykilosort');

    % define save location for quality metrics
    savePath = fullfile(ephys_path, 'qMetrics');
    decompressDataLocal = 'D:\AP006\2023-06-13\'; % where to save raw decompressed ephys data

    %% run/load bombcell
    %% - load data
    [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
        pcFeatureIdx, channelPositions] = bc_loadEphysData(ephys_ks_path);

    %% - detect whether data is compressed, decompress locally if necessary
    rawFile = bc_manageDataCompression(ephys_raw_dir, decompressDataLocal);

    %% - which quality metric parameters to extract and thresholds
    param = bc_qualityParamValues(ephys_meta_dir, rawFile);
    % param = bc_qualityParamValuesForUnitMatch(ephysMetaDir, rawFile) % Run
    % this if you want to use UnitMatch after

    param.nChannels = 384;
    param.nSyncChannels = 0;

    %% - compute quality metrics
    qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

    % I added this manually otherwise it wouldn't extract waveforms and
    % nothing worked
    param.extractRaw = 1;

    if qMetricsExist == 0 || rerun
        [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
            templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,channelPositions, savePath);

        % load raw waveforms
        rawWaveforms.average = readNPY([fullfile(savePath, 'templates._bc_rawWaveforms.npy')]);
        rawWaveforms.peakChan = readNPY([fullfile(savePath, 'templates._bc_rawWaveformPeakChannels.npy')]);

        ephys.extra_bc_q_metrics;

        save(fullfile(savePath, 'am_bc_unit_type.mat'), 'unitType', 'am_bc_units','potential_noise', 'extra_qmetric', '-v7.3')
    else
        [param, qMetric] = bc_loadSavedMetrics(savePath);
%         unitType = bc_getQualityUnitType(param, qMetric);
        load(fullfile(savePath, 'am_bc_unit_type.mat'));
    end

    %% view units + quality metrics in GUI
    % % load data for GUI
%     loadRawTraces = 1;
%     bc_loadMetricsForGUI;
% 
%     % GUI guide:
%     % left/right arrow: toggle between units
%     % g : go to next good unit
%     % m : go to next multi-unit
%     % n : go to next noise unit
%     % up/down arrow: toggle between time chunks in the raw data
%     % u: brings up a input dialog to enter the unit you want to go to
%     unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
%         param, probeLocation, unitType, loadRawTraces);

end
