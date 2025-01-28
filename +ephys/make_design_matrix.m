function design_matrix = make_design_matrix(regressors, t_shifts, time)

design_matrix = [];
for reg_idx=1:length(regressors)
    regressor_times = regressors{reg_idx};
    regressor_t_shift = t_shifts{reg_idx};

    regressor = zeros(size(time));
    timestamps_idx = arrayfun(@(X) find(time < X, 1, 'last'), regressor_times);
    regressor(timestamps_idx) = 1;

    % get idx for stim time shift
    timestep = mean(diff(time));
    t_shift = regressor_t_shift(1):timestep:regressor_t_shift(2);
    t_shift_idx = double(int32(t_shift * 1/timestep));

    % calculate design matrix for stim regressor
    regressor_matrix = lagmatrix(regressor,t_shift_idx);
    regressor_matrix(isnan(regressor_matrix)) = 0;

    design_matrix = horzcat(design_matrix, regressor_matrix);
end