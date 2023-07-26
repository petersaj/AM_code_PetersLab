%% ~~ Example bombcell pipeline ~~
% Adjust the paths in the 'set paths' section and the parameters in bc_qualityParamValues
% This pipeline will (1) load your ephys data, (2) decompress your raw data if it
% is in .cbin format (3) run bombcell on your data and save the output and
% finally (4) bring up summary plots and a GUI to flip through classifyed
% cells.
% The first time, this pipeline will be significantly slower (10-20' more)
% than after because it extracts raw waveforms. Subsequent times these
% pre-extracted waveforms are simply loaded in.
% We recomment running this pipeline on a few datasets and deciding on
% quality metric thresholds depending on the histogram and GUI. 


%% set paths - EDIT THESE 
% ephysKilosortPath = 'P:\Data\AP009\2023-06-27\ephys\pykilosort';% path to your kilosort output files 
% ephysRawDir = dir('D:\AP006\2023-06-13\ephys\probe_1\experiment1\recording1\continuous\Neuropix-3a-100.Neuropix-3a-AP\*.*dat'); % path to yourraw .bin or .dat data
% ephysMetaDir = dir('D:\AP006\2023-06-13\ephys\probe_1\experiment1\recording1\*.*oebin'); % path to your meta / oebin file
% saveLocation = 'D:\AP006\2023-06-13\'; % where you want to save the quality metrics 
% savePath = fullfile(saveLocation, 'qMetrics'); 
% decompressDataLocal = 'D:\AP006\2023-06-13\'; % where to save raw decompressed ephys data 

%% load data 
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc_loadEphysData(ephys_ks_path);

%% detect whether data is compressed, decompress locally if necessary
rawFile = bc_manageDataCompression(ephys_raw_dir, decompressDataLocal);

%% which quality metric parameters to extract and thresholds 
param = bc_qualityParamValues(ephys_meta_dir, rawFile); 
% param = bc_qualityParamValuesForUnitMatch(ephysMetaDir, rawFile) % Run
% this if you want to use UnitMatch after

param.nChannels = 384;
param.nSyncChannels = 0;

%% compute quality metrics 
rerun = 0;
qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

% I added this manually otherwise it wouldn't extract waveforms and
% nothing worked
param.extractRaw = 1;

if qMetricsExist == 0 || rerun
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,channelPositions, savePath);
else
    [param, qMetric] = bc_loadSavedMetrics(savePath); 
    unitType = bc_getQualityUnitType(param, qMetric);
end
