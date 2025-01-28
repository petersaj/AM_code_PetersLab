function [avg_coeff, predict_mua_trace, var_expl] = run_regression(mua_trace, design_matrix, time)

% split into 5
fold_num = 5;
step_for_fold = floor(length(time)/fold_num);
predict_mua_trace = [];
coeff = [];
for fold=1:fold_num
    test_start=(fold - 1) * step_for_fold + 1;
    test_end = test_start + step_for_fold - 1;

    test_idx = false(size(time));
    test_idx(test_start:test_end) = true;
    test_data = mua_trace(test_idx, :);
    test_design_matrix = design_matrix(test_idx, :);

    train_idx = ~test_idx;
    train_data = mua_trace(train_idx,:);
    train_design_matrix = design_matrix(train_idx, :);
% 
%     tic
%     regress_coeff = mvregress(train_design_matrix, train_data);
%     time_for_regression = toc
%     
%     fold_predict_mua_trace = test_design_matrix * regress_coeff;
%     fold_coeff = regress_coeff';

    % using regress
    regress_coeff = regress(train_data, train_design_matrix);
    fold_predict_mua_trace = test_design_matrix * regress_coeff;
    fold_coeff = regress_coeff';
    %     rmse_regress = sqrt(mean((test_data - fitlm_fold_predict_mua_trace).^2));

%     test_coeff = inv(design_matrix*design_matrix')*design_matrix'*test_data

    % OLD %%%%%%%%%%%%%%%
%     % fit linear model
%     lm = fitlm(train_design_matrix, train_data, 'linear', 'Intercept', false);
%     fold_predict_mua_trace = predict(lm, test_design_matrix);
%     %     rmse_fitlm = sqrt(mean((test_data - regress_fold_predict_mua_trace).^2));
% 
%     %      fold_coeff = fold_predict_mua_trace \ test_design_matrix;
% 
%     fold_coeff = lm.Coefficients.Estimate';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % add to big trace
    predict_mua_trace = [predict_mua_trace; fold_predict_mua_trace];
    coeff = [coeff; fold_coeff];
end
avg_coeff = mean(coeff, 1);

% get ssms
var_expl = 1 - (var(mua_trace - predict_mua_trace)) / (var(mua_trace));

%     figure;
%     plot(mua_trace);
%     hold on;
%     plot(predict_mua_trace)
%
%
%     figure;
%     plot(avg_coeff(1:50))
%     hold on;
%     plot(avg_coeff(51:200))
%     hold on;
%     plot(avg_coeff(201:350))
%     hold on;
%     plot(avg_coeff(351:400))
%     legend({'Stim', 'Stim Move', 'ITI Move', 'Reward'})
%     title([rec_day ' Kernels for mua trace str 1'])

