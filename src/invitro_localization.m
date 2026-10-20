% invitro_localization.m
% This script processes 16 microphone audio signals in an air-only medium
% and analyzes the crackle location using DAS, MVDR, and MUSIC algorithms.

% Add the path to mic_preprocess.m (if needed)
% addpath('<path_to_mic_preprocess_folder>');
clear; close all;
% Define a list of folder paths along with optional start and end times
% and actual known source positions.
folder_data = {
    struct('path', '..\data\invitro_02_18_select', 'start', [7.61], 'end', [7.63], 'actual_pos', [-0.05, 0.09, 0.06]); %lower left
    %struct('path', '..\data\invitro_02_18_select', 'start', [12.61], 'end', [12.63], 'actual_pos', [-0.05, 0.02, 0.06]); %lower right
    struct('path', '..\data\invitro_02_18_select', 'start', [16.62], 'end', [16.64], 'actual_pos', [0.02, 0.03, 0.06]); %upper right
    %struct('path', '..\data\invitro_02_18_select', 'start', [22.58], 'end', [22.60], 'actual_pos', [0.03, 0.09, 0.06]); %upper left
    %struct('path', '..\data\2024_11_29_06_21_39_above_mic12_0', 'start', [], 'end', [], 'actual_pos', [0, -0.08, 0.22])
};

%% Step 1: Preprocess Data
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

    if isempty(data.start) || isempty(data.end)
        [cropped_audio_data, common_window_time_ms, latencies_ms] = invitro_mic_preprocess(folder_path);
    else
        [cropped_audio_data, common_window_time_ms, latencies_ms] = invitro_mic_preprocess(folder_path, data.start, data.end);
    end
%% Step 2: Formulate Signal Matrix
%disp('Formulating signal matrix...');
%num_mics = length(cropped_audio_data); % Should be 16
num_mics = 5; %5
signal_matrix = zeros(length(cropped_audio_data{1}), num_mics);

for mic = 1:num_mics
    signal_matrix(:, mic) = cropped_audio_data{mic};
end
    
%% Step 3: Precompute Delays
%disp('Precomputing delays...');
grid_size = 50; % Resolution of the scanning grid
%radius = 0.08; % 8 cm radius around the origin (scaled for air-only setup)
Fs = 64000; % Sampling rate in Hz

% Define scanning grid
x_scan = linspace(-0.12, 0.05, grid_size);
y_scan = linspace(-0.01, 0.11, grid_size);
z_scan = linspace(0.06, 0.06, grid_size);
[X_scan, Y_scan, Z_scan] = ndgrid(x_scan, y_scan, z_scan);
%[X_scan, Y_scan] = ndgrid(x_scan, y_scan);
% Microphone positions in a circle around the origin
%theta = linspace(0, 2*pi, num_mics + 1); theta(end) = [];
%mic_positions = [radius * cos(theta); radius * sin(theta); zeros(1, num_mics)];

mic_positions= [
    -7    -1   0
    -4     8   0
    -1     2   0
     2     5   0
     5    11   0]'/100;
%{
mic_positions= [
    -7    -1   0
     2     5   0
     5    11   0]'/100
%}

% Speed of sound in air
c = 27; % m/s


%% Step 4: SRP-DAS Beamforming 
%disp('Running SRP-DAS...');
das_beamformer = zeros(grid_size, grid_size, grid_size);
for i = 1:numel(X_scan(:))
    % Extract test point
    test_position = [X_scan(i), Y_scan(i), Z_scan(i)];
    combined_signal = zeros(size(signal_matrix, 1), 1);
    
    % Sum contributions from all microphones
    for mic = 1:num_mics
        distance = norm(mic_positions(:, mic) - test_position');
        %delay = round(distance * Fs / c + cal(mic) * Fs/1000); % cal
         delay = round(distance * Fs / c); % cal
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

%% Step 5: Display Results
disp('Beamforming Results:');
disp(['Actual Position: ', num2str(actual_position)]);
disp(['DAS Predicted Position: ', num2str(predicted_position_das'),'   Localization Error: ', num2str(localization_error_das), ' meters']);
%disp(['MVDR Predicted Position: ', num2str(predicted_position_mvdr'), '   Localization Error: ', num2str(localization_error_mvdr), ' meters']);
%disp(['MUSIC Predicted Position: ', num2str(predicted_position_music'),'   Localization Error: ', num2str(localization_error_music), ' meters']);

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
xlabel('X (m)'); xlim([-0.12 0.12]);
ylabel('Y (m)'); ylim([-0.12 0.12]);
zlabel('Z (m)');
grid on;
hold off;

pause(2); % Pause before processing the next folder
end