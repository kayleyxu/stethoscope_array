clear; close all;

%output_filename = '../data/outputs/error_statistics_large_dataset.xlsx';
output_filename = '../data/outputs/error_statistics_small_dataset.xlsx';
%% Load Simulation Data
%load('simulation_results_ml_all_large_dataset/results_for_ml.mat');
load('simulation_results_ml_all_small_dataset/results_for_ml.mat');
data_table = results_for_ml; 

%% Define Dataset Portions
num_features = 6 + 6; % SNR, frequency, c, SRP predicted position (3), and 6 mic delays
num_mics =6;
mic_signal_length = size(data_table, 2) - (num_features + 3 + num_mics * 3); % Extract raw mic signal size
% results_for_ml = [results_for_ml; snr, freq, c, source_position, predicted_position_srp, estimated_delays', reshape(position_sets{1}, 1, []), raw_mic_signals];

%% Extract Features for ML-Only and BF-Assisted ML
X_ml_only = [data_table(:, 4:6), data_table(:, end-mic_signal_length+1:end)]; % ML-Only: Raw mic signals + environmental info
X_bf_assisted = [data_table(:, 4:6), data_table(:, 7:9), data_table(:, 10:10+num_mics-1)]; % BF-Assisted ML: SRP & TDOA features

% Labels (True source location)
Y = data_table(:, 1:3); % True x, y, z coordinates

%% Split Data: 80% Training, 20% Testing
rng(7); % For reproducibility
% defines a random partition on a data set, 20% test
cv = cvpartition(size(X_ml_only, 1), 'HoldOut', 0.2); 
% after partition, train(cv) returns "logical array" of 1s, 0s to mark
% which index is for training, then X_ml_only(training(cv), :) extracts
% only training entries. Similar for test(cv).
% _ml and _bf can share the same partition and indexing (cv), and same Y
X_train_ml = X_ml_only(training(cv), :);
X_test_ml = X_ml_only(test(cv), :);
X_train_bf = X_bf_assisted(training(cv), :);
X_test_bf = X_bf_assisted(test(cv), :);
Y_train = Y(training(cv), :);
Y_test = Y(test(cv), :);

% assign same corresponding test index to SRP only for fair comparison
srp_only = data_table(:, 7:9); % SRP-DAS only
test_srp = srp_only(test(cv), :);

%% Train Machine Learning Models for ML-Only and BF-Assisted ML
fprintf('Training models...\n');

models = {'LSBoost', 'GPR', 'SVM'};
predictions = struct();
errors = struct();

% Loop through different models' train and test
% ML train and create separate models for x, y, z coordinates
for i = 1:length(models)
    model_name = models{i};
    
    if strcmp(model_name, 'LSBoost')
        train_func = @(X, Y) fitrensemble(X, Y, 'Method', 'LSBoost', 'NumLearningCycles', 100);
    elseif strcmp(model_name, 'GPR')
        train_func = @(X, Y) fitrgp(X, Y);
    elseif strcmp(model_name, 'SVM')
        train_func = @(X, Y) fitrsvm(X, Y);
    end
    
    % Train models for ML-Only
    model_ml_x = train_func(X_train_ml, Y_train(:, 1));
    model_ml_y = train_func(X_train_ml, Y_train(:, 2));
    model_ml_z = train_func(X_train_ml, Y_train(:, 3));
    
    % Train models for BF-Assisted ML
    model_bf_x = train_func(X_train_bf, Y_train(:, 1));
    model_bf_y = train_func(X_train_bf, Y_train(:, 2));
    model_bf_z = train_func(X_train_bf, Y_train(:, 3));
    
    % Predict on Test Data
    predictions.(model_name).ml = [predict(model_ml_x, X_test_ml), predict(model_ml_y, X_test_ml), predict(model_ml_z, X_test_ml)];
    predictions.(model_name).bf = [predict(model_bf_x, X_test_bf), predict(model_bf_y, X_test_bf), predict(model_bf_z, X_test_bf)];
    
    % Compute Errors
    errors.(model_name).ml = vecnorm(predictions.(model_name).ml - Y_test, 2, 2);
    errors.(model_name).bf = vecnorm(predictions.(model_name).bf - Y_test, 2, 2);
end

%% Compute Errors for SRP-DAS Only
error_srp = vecnorm(test_srp - Y_test, 2, 2);

%{
%% Compare Results
fprintf('Mean SRP-DAS Error: %.4f\n', mean(error_srp));
for i = 1:length(models)
    model_name = models{i};
    fprintf('Mean %s ML-Only Error: %.4f\n', model_name, mean(errors.(model_name).ml));
    fprintf('Mean %s BF-Assisted Error: %.4f\n', model_name, mean(errors.(model_name).bf));
end
%}

%% Statistical Analysis
errors_table = array2table([
    mean(error_srp), median(error_srp), sqrt(mean(error_srp.^2)), var(error_srp);
    mean(errors.LSBoost.ml), median(errors.LSBoost.ml), sqrt(mean(errors.LSBoost.ml.^2)), var(errors.LSBoost.ml);
    mean(errors.GPR.ml), median(errors.GPR.ml), sqrt(mean(errors.GPR.ml.^2)), var(errors.GPR.ml);
    mean(errors.SVM.ml), median(errors.SVM.ml), sqrt(mean(errors.SVM.ml.^2)), var(errors.SVM.ml);
    mean(errors.LSBoost.bf), median(errors.LSBoost.bf), sqrt(mean(errors.LSBoost.bf.^2)), var(errors.LSBoost.bf);
    mean(errors.GPR.bf), median(errors.GPR.bf), sqrt(mean(errors.GPR.bf.^2)), var(errors.GPR.bf);
    mean(errors.SVM.bf), median(errors.SVM.bf), sqrt(mean(errors.SVM.bf.^2)), var(errors.SVM.bf)
], 'VariableNames', {'Mean', 'Median', 'RMSE', 'Variance'}, ...
   'RowNames', {'SRP-DAS', 'LSBoost-ML', 'GPR-ML', 'SVM-ML', 'LSBoost-BF', 'GPR-BF', 'SVM-BF'});

disp('Error Comparison Table:');
disp(errors_table);

% Export table to an Excel file
writetable(errors_table, output_filename, 'WriteRowNames', true);

%% Prepare Data for Violin Plot (No Padding Needed)
violindata = [error_srp, errors.LSBoost.ml, errors.GPR.ml, errors.SVM.ml, ...
              errors.LSBoost.bf, errors.GPR.bf, errors.SVM.bf];

%% Generate Visualizations
figure;

% Boxplot with Custom Colors
subplot(1,2,1);
%colors = {[0 0 1], [0 0 1], [0 0 1], [0 1 0], [0 1 0], [0 1 0],[1 0 0]}; % SRP-DAS (blue), ML-Only (green), BF-ML (red)
colors = {'g','g','g','m','m','m','b'} % SRP-DAS (blue), ML-Only, BF-ML (green)
box = boxplot(violindata, 'Labels', {'SRP-DAS', 'LSBoost-ML', 'GPR-ML', 'SVM-ML', 'LSBoost-BF', 'GPR-BF', 'SVM-BF'});
title('Error Distribution (Boxplot)');
ylabel('Error');
grid on;
ylim([0, max(violindata(:)) * 1.1]);

% Apply Different Colors for Different Groups
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{j}, 'FaceAlpha', .5);
end

% Violin Plot with Proper X-axis Labels and Y-limit Fix
subplot(1,2,2);
v = violinplot(violindata);
xticklabels({'SRP-DAS', 'LSBoost-ML', 'GPR-ML', 'SVM-ML', 'LSBoost-BF', 'GPR-BF', 'SVM-BF'});
title('Error Distribution (Violin Plot)');
ylabel('Error');
grid on;
ylim([0, max(violindata(:)) * 1.1]); % Ensure Y-axis is strictly positive

sgtitle('Error Comparison Across Methods');
fig = gcf;  % Get the current figure handle
savefig(fig, '../data/outputs/error_distribution_large_dataset.fig');
disp('Figure saved as error_distribution.fig');
