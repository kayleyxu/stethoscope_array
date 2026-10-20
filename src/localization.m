% Processes 16 microphone audio signals in an air-only medium
% and analyzes the sound (crackle) location using SRP-DAS algorithm

% Add the path to mic_preprocess.m (if needed)
% addpath('<path_to_mic_preprocess_folder>');
clear; close all;
% Define the file path
%folder_path = '..\data\2024_11_26_23_51_25_single_crackle';
%folder_path = '..\data\2024_11_29_05_55_26_center_0'; %2,3 good
%folder_path = '..\data\2024_11_29_06_21_39_center_1'; %6.5, 7.6
%folder_path = '..\data\2024_11_29_06_21_39_above_mic8_0'; %2.5,3.5 put in range good
%folder_path = '..\data\2024_11_29_06_21_39_above_mic8_1'; %4,5 
%folder_path = '..\data\2024_11_29_06_21_39_above_mic0_0';
%folder_path = '..\data\2024_11_29_06_21_39_above_mic4_0';
%folder_path = '..\data\2024_11_29_06_21_39_above_mic12_0';
% Define a list of folder paths along with optional start and end times
% and actual known source positions.
folder_data = {
    struct('path', '..\data\2024_11_26_23_51_25_single_crackle', 'start', [], 'end', [], 'actual_pos', [-0.08, 0.0, 0.22]);
    struct('path', '..\data\2024_11_29_05_55_26_center_0', 'start', 2, 'end', 3, 'actual_pos', [0.0, 0.0, 0.22]);
    %struct('path', '..\data\2024_11_29_06_21_39_center_1', 'start', 6.5, 'end', 7.6, 'actual_pos', [0.0, 0.0, 0.22]);
    struct('path', '..\data\2024_11_29_06_21_39_above_mic0_0', 'start', [], 'end', [], 'actual_pos', [0.08, 0.0, 0.22]);
    struct('path', '..\data\2024_11_29_06_21_39_above_mic4_0', 'start', [], 'end', [], 'actual_pos', [0, 0.08, 0.22]);
    struct('path', '..\data\2024_11_29_06_21_39_above_mic8_0', 'start', 2.5, 'end', 3.5, 'actual_pos', [-0.08, 0, 0.22]);
    %struct('path', '..\data\2024_11_29_06_21_39_above_mic8_1', 'start', 4, 'end', 5, 'actual_pos', [-0.08, 0, 0.22]);
    struct('path', '..\data\2024_11_29_06_21_39_above_mic12_0', 'start', [], 'end', [], 'actual_pos', [0, -0.08, 0.22])
};

%% Preprocess Data
%disp(['Preprocessing data from ', folder_path]);
%[cropped_audio_data, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path);
% Display the start time of the common window
%disp(['Using common 20ms window starting at: ', num2str(round(common_window_time_ms, 2)), ' ms']);
% Loop through each folder and process the data
for i = 1:length(folder_data)
    data = folder_data{i};
    folder_path = data.path; 
    actual_position = data.actual_pos;

    disp(['Processing folder: ', folder_path]);

    %% Preprocess Data (with optional start and end times)
    if isempty(data.start) || isempty(data.end)
        [cropped_audio_data, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path);
    else
        [cropped_audio_data, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path, data.start, data.end);
    end
%% Step 2: Formulate Signal Matrix
%disp('Formulating signal matrix...');
num_mics = length(cropped_audio_data); % Should be 16
signal_matrix = zeros(length(cropped_audio_data{1}), num_mics);

for mic = 1:num_mics
    signal_matrix(:, mic) = cropped_audio_data{mic};
end

%% Precompute Delays
%disp('Precomputing delays...');
grid_size = 50; % Resolution of the scanning grid
radius = 0.08; % 8 cm radius around the origin (scaled for air-only setup)
Fs = 48000; % Sampling rate in Hz

% Define scanning grid
x_scan = linspace(-radius, radius, grid_size);
y_scan = linspace(-radius, radius, grid_size);
z_scan = linspace(0.22, 0.22, grid_size);
[X_scan, Y_scan, Z_scan] = ndgrid(x_scan, y_scan, z_scan);
%[X_scan, Y_scan] = ndgrid(x_scan, y_scan);
% Microphone positions in a circle around the origin
theta = linspace(0, 2*pi, num_mics + 1); theta(end) = [];
mic_positions = [radius * cos(theta); radius * sin(theta); zeros(1, num_mics)];

% Speed of sound in air
c_air = 343; % m/s

%% SRP-DAS Beamforming 
%disp('Running SRP-DAS...');
das_beamformer = zeros(grid_size, grid_size, grid_size);
for i = 1:numel(X_scan(:))
    % Extract test point
    test_position = [X_scan(i), Y_scan(i), Z_scan(i)];
    combined_signal = zeros(size(signal_matrix, 1), 1);
    
    % Sum contributions from all microphones
    for mic = 1:num_mics
        distance = norm(mic_positions(:, mic) - test_position');
        %delay = round(distance * Fs / c_air + cal(mic) * Fs/1000); % cal
         delay = round(distance * Fs / c_air); % cal
        % Time-shift and sum signals
        shifted_signal = circshift(signal_matrix(:, mic), -delay);
        combined_signal = combined_signal + shifted_signal;
    end

    % Compute energy for the current test point
    das_beamformer(i) = sum(combined_signal.^2);
end

% Find maximum point in DAS beamformer
[~, max_idx] = max(das_beamformer(:));
[grid_x, grid_y, grid_z] = ind2sub(size(das_beamformer), max_idx);
predicted_position_das = [x_scan(grid_x); y_scan(grid_y); z_scan(grid_z)];
localization_error_das = norm(predicted_position_das' - actual_position);

%% Display Results
disp('Beamforming Results:');
disp(['Actual Position: ', num2str(actual_position)]);
disp(['DAS Predicted Position: ', num2str(predicted_position_das'),'   Localization Error: ', num2str(localization_error_das), ' meters']);
%disp(['MVDR Predicted Position: ', num2str(predicted_position_mvdr'), '   Localization Error: ', num2str(localization_error_mvdr), ' meters']);
%disp(['MUSIC Predicted Position: ', num2str(predicted_position_music'),'   Localization Error: ', num2str(localization_error_music), ' meters']);

%% Visualize Predicted Positions
figure;
scatter3(mic_positions(1, :), mic_positions(2, :), mic_positions(3, :), 'g', 'filled');
hold on;

% Label mic indices
for mic = 1:num_mics
    text(mic_positions(1, mic), mic_positions(2, mic), mic_positions(3, mic), ...
        ['Mic ', num2str(mic-1)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

scatter3(actual_position(1), actual_position(2), actual_position(3), 100,'rd', 'filled');
scatter3(predicted_position_das(1), predicted_position_das(2), predicted_position_das(3), 'b', 'filled');
%scatter3(predicted_position_mvdr(1), predicted_position_mvdr(2), predicted_position_mvdr(3), 'c', 'filled');
%scatter3(predicted_position_music(1), predicted_position_music(2), predicted_position_music(3), 'm', 'filled');

%legend('Microphones', 'Actual', 'DAS', 'MVDR', 'MUSIC');
legend('Microphones', 'Actual', 'DAS Predicted');
%legend('Microphones', 'MVDR', 'MUSIC');
title('Predicted Sound Source Positions');
xlabel('X (m)'); xlim([-0.1 0.1]);
ylabel('Y (m)'); ylim([-0.1 0.1]);
zlabel('Z (m)');
grid on;
hold off;

pause(1); % Pause before processing the next folder
end