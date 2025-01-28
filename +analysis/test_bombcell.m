savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\new_qMetrics';

param = parquetread([fullfile(savePath, '_bc_parameters._bc_qMetrics.parquet')]);

param.minWidthFirstPeak %= 4; % in samples, minimum width of waveform's first peak
param.minMainPeakToTroughRatio % = 10; % min peak to trough ratio
param.minWidthMainTrough %= 5; % in samples, minimum width of waveform's main trough

[param, qMetric] = bc.load.loadSavedMetrics(savePath); 

% doesn't exist qMetric.minWidthMainTrough

% tried just changing the paramater to see if it does anything - NOPE

isNonSomatic = (abs(qMetric.mainTrough_size ./ qMetric.mainPeak_before_size) < param.minMainPeakToTroughRatio &...
        qMetric.mainPeak_before_width < param.minWidthFirstPeak &...
        qMetric.mainTrough_width < param.minWidthMainTrough &...
        abs(qMetric.mainPeak_before_size./qMetric.mainPeak_after_size) > param.firstPeakRatio) | ...
        abs(max([qMetric.mainPeak_before_size, qMetric.mainPeak_after_size], [], 2)./ qMetric.mainTrough_size) > param.minTroughToPeakRatio;

abs(max([qMetric.mainPeak_before_size(unit_idx), qMetric.mainPeak_after_size(unit_idx)], [], 2)./ qMetric.mainTrough_size(unit_idx))
qMetric.mainTrough_width(unit_idx)


unit_idx = 752; % axon that pases as mua
unit_idx = 764; % axon that pases as single unit
unit_idx = 629; % mua
unit_idx = 591; % mua
unit_idx = 572; % mua

% for mainThrough_width: 629, 591 and 572 have 0, (width around 6) 
% both of the axons have 1 (width 3/4)
% ---> change the logic of the statement here maybe?? 
(abs(qMetric.mainTrough_size(unit_idx) ./ qMetric.mainPeak_before_size(unit_idx)) < param.minMainPeakToTroughRatio &...
        qMetric.mainPeak_before_width(unit_idx) < param.minWidthFirstPeak &...
        qMetric.mainTrough_width(unit_idx) < param.minWidthMainTrough &...
        abs(qMetric.mainPeak_before_size(unit_idx)./qMetric.mainPeak_after_size(unit_idx)) > param.firstPeakRatio) | ...
        abs(max([qMetric.mainPeak_before_size(unit_idx), qMetric.mainPeak_after_size(unit_idx)], [], 2)./ qMetric.mainTrough_size(unit_idx)) > param.minTroughToPeakRatio

% like this
new_isNonSomatic = (abs(qMetric.mainTrough_size ./ qMetric.mainPeak_before_size) < param.minMainPeakToTroughRatio &...
        qMetric.mainPeak_before_width < param.minWidthFirstPeak &...
        abs(qMetric.mainPeak_before_size./qMetric.mainPeak_after_size) > param.firstPeakRatio) | ...
        abs(max([qMetric.mainPeak_before_size, qMetric.mainPeak_after_size], [], 2)./ qMetric.mainTrough_size) > param.minTroughToPeakRatio | ...
        qMetric.mainTrough_width < param.minWidthMainTrough;

axon_idx_org = isNonSomatic & ~new_isNonSomatic;
axon_idx_new = ~isNonSomatic & new_isNonSomatic;

% nope, this gives me almost no units :( 

narrow_axons = qMetric.mainTrough_width < 4 & ...
    qMetric.mainPeak_before_width < 3 ;


find(unitType == 5)

%% plot
% - param.minMainPeakToTroughRatio and
% - param.firstPeakRatio

figure;
scatter(abs(qMetric.mainTrough_size ./ qMetric.mainPeak_before_size), ...
    abs(qMetric.mainPeak_before_size./qMetric.mainPeak_after_size)) 
xline(param.firstPeakRatio)
yline(param.minMainPeakToTroughRatio)
xlabel('minMainPeakToTroughRatio')
ylabel('firstPeakRatio')

% through width vs peak before width
figure;
scatter(qMetric.mainTrough_width, ...
     qMetric.mainPeak_before_width) 
xline(4)
yline(4)
xlabel('mainTrough width')
ylabel('mainPeak before width')

% 
figure;
scatter(qMetric.mainTrough_width, ...
    abs(qMetric.mainPeak_before_size./qMetric.mainPeak_after_size)) 
xline(4)
yline(param.firstPeakRatio)
xlabel('mainTrough width')
ylabel('firstPeakRatio')

%% new metrics: width 4 peak 4

savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\width_4_and_peak_4_qMetrics';

param = parquetread([fullfile(savePath, '_bc_parameters._bc_qMetrics.parquet')]);

% param.minWidthFirstPeak %= 4; % in samples, minimum width of waveform's first peak
% param.minMainPeakToTroughRatio % = 5; % min peak to trough ratio
% param.minWidthMainTrough %= 5; % in samples, minimum width of waveform's main trough
% param.mainPeak_before_width 

[param, qMetric] = bc.load.loadSavedMetrics(savePath); 

unit_idx = 752; % axon 
unit_idx = 764; % axon 
unit_idx = 629; % mua
unit_idx = 591; % mua
unit_idx = 572; % mua

unitType(unit_idx)

%% for old metrics
savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\qMetrics';

param = parquetread([fullfile(savePath, '_bc_parameters._bc_qMetrics.parquet')]);

[param, qMetric] = bc.load.loadSavedMetrics(savePath); 

animal = 'AM014';
rec_day = '2024-01-23';
workflow = 'passive';
load_parts.ephys = true;
ap.load_recording

trial_stim_values = vertcat(trial_events.values.TrialStimX);
trial_stim_values = trial_stim_values(1:length(stimOn_times));

ap.cellraster(stimOn_times,trial_stim_values);


figure;
scatter3(qMetric.mainTrough_width(good_templates), ...
    abs(qMetric.mainPeak_after_size(good_templates)./qMetric.mainPeak_before_size(good_templates)), ...
    qMetric.mainPeak_before_width(good_templates)) 



% xline(4)
% yline(param.firstPeakRatio)
xlabel('mainTrough width')
ylabel('firstPeakRatio')
zlabel('first peak width')


%% compare
load("P:\Data\AM014\2024-01-23\ephys\kilosort4\qMetrics\template_qc_labels.mat")
rerun_labels = template_qc_labels;

load("P:\Data\AM014\2024-01-23\ephys\kilosort4\new_qMetrics\template_qc_labels.mat")
new_labels = template_qc_labels;

setdiff(rerun_labels, new_labels)

diff_noise_idx = find(strcmpi(rerun_labels, new_labels)==0);
rerun_labels(diff_noise_idx)
new_labels(diff_noise_idx)

% rerun
rerun_savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\qMetrics';
rerun_param = parquetread([fullfile(rerun_savePath, '_bc_parameters._bc_qMetrics.parquet')]);
[rerun_param, rerun_qMetric] = bc.load.loadSavedMetrics(rerun_savePath); 
% new
new_savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\new_qMetrics';
new_param = parquetread([fullfile(new_savePath, '_bc_parameters._bc_qMetrics.parquet')]);
[new_param, new_qMetric] = bc.load.loadSavedMetrics(new_savePath); 

% compare
param_same = isequaln(rerun_param, new_param);
qMetric_same = isequaln(new_qMetric, rerun_qMetric);


% Step 1: Check dimensions - same
size1 = size(rerun_qMetric);
size2 = size(new_qMetric);
if ~isequal(size1, size2)
    disp('The tables have different dimensions:');
    disp(['rerun_qMetric size: ', mat2str(size1)]);
    disp(['new_qMetric size: ', mat2str(size2)]);
end

% Step 2: Compare variable names - same
varNames1 = rerun_qMetric.Properties.VariableNames;
varNames2 = new_qMetric.Properties.VariableNames;

if ~isequal(varNames1, varNames2)
    disp('The tables have different variable names:');
    disp('rerun_qMetric variable names:');
    disp(varNames1);
    disp('new_qMetric variable names:');
    disp(varNames2);
end

% Step 3: Compare data column-wise - all the same
commonVars = intersect(varNames1, varNames2); % Variables present in both tables
for i = 1:length(commonVars)
    varName = commonVars{i};
    col1 = rerun_qMetric.(varName);
    col2 = new_qMetric.(varName);
    
    if isnumeric(col1) || islogical(col1)
        % Handle NaN: Treat NaNs in the same position as equal
        nanMask = isnan(col1) & isnan(col2);
        col1(nanMask) = 0; % Replace NaN with a dummy value
        col2(nanMask) = 0; % Replace NaN with the same dummy value
        differences = col1 ~= col2;
    elseif iscell(col1) || isstring(col1) || ischar(col1)
        % Handle string or cell arrays
        differences = ~strcmp(col1, col2);
    else
        disp(['Column "', varName, '" data type is not directly comparable.']);
        continue;
    end

    % Display differences
    if any(differences)
        disp(['Column "', varName, '" differs between the tables.']);
%         diffIndices = find(differences);
%         disp(table(diffIndices, col1(differences), col2(differences), ...
%             'VariableNames', {'RowIndex', 'rerun_qMetric', 'new_qMetric'}));
    end
end

% Step 4: Check for missing columns
missingInRerun = setdiff(varNames2, varNames1);
missingInNew = setdiff(varNames1, varNames2);

if ~isempty(missingInRerun)
    disp('Columns in new_qMetric but not in rerun_qMetric:');
    disp(missingInRerun);
end

if ~isempty(missingInNew)
    disp('Columns in rerun_qMetric but not in new_qMetric:');
    disp(missingInNew);
end

% check how many in column
% Define the column name to check
columnName = 'rawAmplitude';

% Check if the column exists in both tables
if ismember(columnName, rerun_qMetric.Properties.VariableNames) && ...
   ismember(columnName, new_qMetric.Properties.VariableNames)
    
    % Extract the columns from each table
    col1 = rerun_qMetric.(columnName);
    col2 = new_qMetric.(columnName);
    
    % Handle NaN equivalence and find differences
    if isnumeric(col1) || islogical(col1)
        % Treat NaN values as equal
        nanMask = isnan(col1) & isnan(col2);
        col1(nanMask) = 0; % Replace NaN with a dummy value
        col2(nanMask) = 0;
        differences = col1 ~= col2;
    elseif iscell(col1) || isstring(col1) || ischar(col1)
        % For string or cell arrays, use strcmp to find differences
        differences = ~strcmp(col1, col2);
    else
        error('Unsupported data type in the specified column.');
    end
    
    % Count the number of differing rows
    numDifferences = sum(differences);
    disp(['Number of differing rows in column "', columnName, '": ', num2str(numDifferences)]);
    
    % Optionally display the indices of differing rows
%     if numDifferences > 0
%         disp('Indices of differing rows:');
%         disp(find(differences));
%     end
else
    disp(['Column "', columnName, '" does not exist in one or both tables.']);
end

% HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rerun_rawWaveformsFull = readNPY(fullfile(rerun_savePath, 'templates._bc_rawWaveforms_kilosort_format.npy'));
new_rawWaveformsFull = readNPY(fullfile(new_savePath, 'templates._bc_rawWaveforms_kilosort_format.npy'));

rerun_rawWaveformsPeakChan = readNPY(fullfile(rerun_savePath, 'templates._bc_rawWaveformPeakChannels_kilosort_format.npy'));
new_rawWaveformsPeakChan = readNPY(fullfile(new_savePath, 'templates._bc_rawWaveformPeakChannels_kilosort_format.npy'));

diff_idx = find(differences);

figure;
tiledlayout('flow')
for this_idx=1:length(diff_idx(1:10))
    i = diff_idx(this_idx);
    nexttile;
    this_new_rawWaveform = squeeze(new_rawWaveformsFull(i, new_rawWaveformsPeakChan(i), :));
    this_rerun_rawWaveform = squeeze(rerun_rawWaveformsFull(i, rerun_rawWaveformsPeakChan(i), :));
    this_same_rerun_rawWaveform = squeeze(rerun_rawWaveformsFull(i, new_rawWaveformsPeakChan(i), :));
    plot(this_new_rawWaveform)
    hold on;
    plot(this_rerun_rawWaveform)
    hold on;
    plot(this_same_rerun_rawWaveform)
    title(['NEW: ' num2str(new_rawWaveformsPeakChan(i)) ...
        '; RERUN: ' num2str(rerun_rawWaveformsPeakChan(i)) ...
        '; SAME: ' num2str(new_rawWaveformsPeakChan(i))])
end

% load old
load("P:\Data\AM014\2024-01-23\ephys\kilosort4\old_qMetrics\template_qc_labels.mat")
old_labels = template_qc_labels;

old_savePath = 'P:\Data\AM014\2024-01-23\ephys\kilosort4\old_qMetrics';
old_param = parquetread([fullfile(old_savePath, '_bc_parameters._bc_qMetrics.parquet')]);
[old_param, old_qMetric] = bc.load.loadSavedMetrics(old_savePath); 

old_rawWaveformsFull = readNPY(fullfile(old_savePath, 'templates._bc_rawWaveforms.npy'));
old_rawWaveformsPeakChan = readNPY(fullfile(old_savePath, 'templates._bc_rawWaveformPeakChannels.npy'));

figure;
tiledlayout('flow')
for this_idx=1:length(diff_idx(1:10))
    i = diff_idx(this_idx);
    nexttile;
    this_new_rawWaveform = squeeze(new_rawWaveformsFull(i, new_rawWaveformsPeakChan(i), :));
    this_rerun_rawWaveform = squeeze(rerun_rawWaveformsFull(i, rerun_rawWaveformsPeakChan(i), :));
    this_old_rawWaveform = squeeze(old_rawWaveformsFull(i, old_rawWaveformsPeakChan(i), :));
    plot(this_new_rawWaveform)
    hold on;
    plot(this_rerun_rawWaveform)
    hold on;
    plot(this_same_rerun_rawWaveform)
    title(['NEW: ' num2str(new_rawWaveformsPeakChan(i)) ...
        '; RERUN: ' num2str(rerun_rawWaveformsPeakChan(i)) ...
        '; OLD: ' num2str(old_rawWaveformsPeakChan(i))])
end

test_rerun_rawWaveformsFull = readNPY(fullfile(rerun_savePath, 'templates._bc_rawWaveforms.npy'));
test_new_rawWaveformsFull = readNPY(fullfile(new_savePath, 'templates._bc_rawWaveforms.npy'));

test_rerun_rawWaveformsPeakChan = readNPY(fullfile(rerun_savePath, 'templates._bc_rawWaveformPeakChannels.npy'));
test_new_rawWaveformsPeakChan = readNPY(fullfile(new_savePath, 'templates._bc_rawWaveformPeakChannels.npy'));

figure;
tiledlayout('flow')
for this_idx=1:length(diff_idx(1:10))
    i = diff_idx(this_idx);
    nexttile;
    this_new_rawWaveform = squeeze(test_new_rawWaveformsFull(i, test_new_rawWaveformsPeakChan(i), :));
%     this_rerun_rawWaveform = squeeze(test_rerun_rawWaveformsFull(i, test_rerun_rawWaveformsPeakChan(i), :));
    this_old_rawWaveform = squeeze(old_rawWaveformsFull(i, test_new_rawWaveformsPeakChan(i), :));
    plot(this_new_rawWaveform)
    hold on;
    plot(this_old_rawWaveform)
    %     hold on;
%     plot(this_rerun_rawWaveform)
    title(['NEW: ' num2str(new_rawWaveformsPeakChan(i)) ...
        '; OLD: ' num2str(old_rawWaveformsPeakChan(i))])
%         '; RERUN: ' num2str(rerun_rawWaveformsPeakChan(i)) ...
        
end

figure;
tiledlayout('flow')
for this_idx=1:length(diff_idx(1:10))
    i = diff_idx(this_idx);
    nexttile;
    this_new_rawWaveform = squeeze(new_rawWaveformsFull(i, new_rawWaveformsPeakChan(i), :));
%     this_rerun_rawWaveform = squeeze(test_rerun_rawWaveformsFull(i, test_rerun_rawWaveformsPeakChan(i), :));
    this_old_rawWaveform = squeeze(old_rawWaveformsFull(i, new_rawWaveformsPeakChan(i), :));
    plot(this_new_rawWaveform)
    hold on;
    plot(this_old_rawWaveform)
    %     hold on;
%     plot(this_rerun_rawWaveform)
    title(['NEW: ' num2str(new_rawWaveformsPeakChan(i)) ...
        '; OLD: ' num2str(old_rawWaveformsPeakChan(i))])
%         '; RERUN: ' num2str(rerun_rawWaveformsPeakChan(i)) ...
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % check data types
% varNames1 = rerun_param.Properties.VariableNames;
% varNames2 = new_param.Properties.VariableNames;
% 
% for i = 1:length(varNames1)
%     varName = varNames1{i};
%     type1 = class(rerun_param.(varName));
%     type2 = class(new_param.(varName));
%     if ~strcmp(type1, type2)
%         disp(['Column "', varName, '" has different data types:']);
%         disp(['rerun_param: ', type1, ', new_param: ', type2]);
%     end
% end
% 
% % Check column order
% if ~isequal(rerun_param.Properties.VariableNames, new_param.Properties.VariableNames)
%     disp('The columns are in a different order.');
% end
% 
% % Check row order (if applicable)
% if ~isequal(rerun_param.Properties.RowNames, new_param.Properties.RowNames)
%     disp('The rows are in a different order.');
% end
% 
% % metadata
% if ~isequal(rerun_param.Properties.Description, new_param.Properties.Description)
%     disp('The tables have different descriptions.');
% end
% 
% if ~isequal(rerun_param.Properties.UserData, new_param.Properties.UserData)
%     disp('The tables have different UserData.');
% end
% 
% % Check custom properties
% customProps1 = fieldnames(rerun_param.Properties.CustomProperties);
% customProps2 = fieldnames(new_param.Properties.CustomProperties);
% 
% if ~isequal(customProps1, customProps2)
%     disp('The tables have different custom properties:');
%     disp('rerun_param custom properties:');
%     disp(customProps1);
%     disp('new_param custom properties:');
%     disp(customProps2);
% end
% 
% % nans
% for i = 1:length(varNames1)
%     varName = varNames1{i};
%     col1 = rerun_param.(varName);
%     col2 = new_param.(varName);
%     if ismissing(col1) ~= ismissing(col2)
%         disp(['Column "', varName, '" has differing missing value patterns.']);
%     end
% end
% 
% % Hidden Differences
% struct1 = table2struct(rerun_param);
% struct2 = table2struct(new_param);
% 
% if ~isequal(struct1, struct2)
%     disp('The tables have subtle differences in their data or structure.');
% end
% 
% disp('Structure of rerun_param:');
% whos rerun_param
% 
% disp('Structure of new_param:');
% whos new_param
% 
% struct1 = table2struct(rerun_param);
% struct2 = table2struct(new_param);
% 
% if ~isequaln(struct1, struct2)
%     disp('Detailed struct comparison shows differences.');
%     % Optionally, use MATLAB's `compare` utility or write a loop to inspect fields
% end


