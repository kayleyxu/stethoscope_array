

clc; clear; close all;

%% Parameters
speed_of_sound = 40; % Speed of sound (m/s)
fs = 64000; % Sampling frequency (Hz)
event_duration = 0.1; % Duration of the event (s)
event_amplitude = 1; % Amplitude of the event
radius = 0.1; % Radius of the circular array in meters
num_sources = 10; % Number of random sound sources per run
num_runs = 100; % Number of runs for averaging
snr_values = [-5];%, 0, 5, 10, 20]; % Noise levels (SNR in dB)

% Manual angles
angles_golomb = deg2rad([0, 17.1429, 68.5714, 240.0000, 274.2857]); % Convert to radians
angles_even = deg2rad([0, 72, 144, 216, 288]); % Convert to radians

% Compute microphone positions for Golomb and evenly spaced arrays
mics_golomb = radius * [cos(angles_golomb)', sin(angles_golomb)'];
mics_even = radius * [cos(angles_even)', sin(angles_even)'];

% Time vector for the event signal
time_vector = 0:1/fs:event_duration;
% Event signal: Gaussian-modulated pulse
fc = 50; % Center frequency (Hz)
sigma = 0.01; % Pulse duration (seconds)
event_signal = event_amplitude * exp(-((time_vector - event_duration/2) / sigma).^2) .* sin(2 * pi * fc * time_vector);

%% Initialize Statistics Storage
avg_errors_golomb = zeros(length(snr_values), 1);
median_errors_golomb = zeros(length(snr_values), 1);
variance_errors_golomb = zeros(length(snr_values), 1);

avg_errors_even = zeros(length(snr_values), 1);
median_errors_even = zeros(length(snr_values), 1);
variance_errors_even = zeros(length(snr_values), 1);

%% Perform Multiple Runs
for snr_idx = 1:length(snr_values)
    snr = snr_values(snr_idx);
    
    errors_golomb = zeros(num_runs, num_sources);
    errors_even = zeros(num_runs, num_sources);
    
    for run = 1:num_runs
        % Generate random source positions
        source_positions = radius * 2 * (rand(num_sources, 2) - 0.5); % Random (x, y) within [-radius, radius]

        for i = 1:num_sources
            % True source position
            source = source_positions(i, :);

            % Simulate signals for Golomb array
            signals_golomb = simulate_propagation(mics_golomb, source, event_signal, fs, speed_of_sound);
            signals_golomb = awgn(signals_golomb, snr, 'measured'); % Add noise
            estimated_golomb = grid_search_localization(mics_golomb, signals_golomb, radius, event_signal, fs, speed_of_sound);
            errors_golomb(run, i) = norm(estimated_golomb - source);

            % Simulate signals for Even array
            signals_even = simulate_propagation(mics_even, source, event_signal, fs, speed_of_sound);
            signals_even = awgn(signals_even, snr, 'measured'); % Add noise
            estimated_even = grid_search_localization(mics_even, signals_even, radius, event_signal, fs, speed_of_sound);
            errors_even(run, i) = norm(estimated_even - source);
        end
    end
    
    % Compute statistics for Golomb array
    all_errors_golomb = errors_golomb(:);
    avg_errors_golomb(snr_idx) = mean(all_errors_golomb);
    median_errors_golomb(snr_idx) = median(all_errors_golomb);
    variance_errors_golomb(snr_idx) = var(all_errors_golomb);

    % Compute statistics for Even array
    all_errors_even = errors_even(:);
    avg_errors_even(snr_idx) = mean(all_errors_even);
    median_errors_even(snr_idx) = median(all_errors_even);
    variance_errors_even(snr_idx) = var(all_errors_even);
end

%% Display Results
disp('Golomb Array Performance:');
disp('SNR | Avg Error | Median Error | Variance');
disp([snr_values', avg_errors_golomb, median_errors_golomb, variance_errors_golomb]);

disp('Evenly Spaced Array Performance:');
disp('SNR | Avg Error | Median Error | Variance');
disp([snr_values', avg_errors_even, median_errors_even, variance_errors_even]);

%% Plot Results
figure;
subplot(1, 3, 1);
plot(snr_values, avg_errors_golomb, '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(snr_values, avg_errors_even, '-rx', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)');
ylabel('Average Error (m)');
title('Average Localization Error');
legend('Golomb Array', 'Evenly Spaced Array');
grid on;

subplot(1, 3, 2);
plot(snr_values, median_errors_golomb, '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(snr_values, median_errors_even, '-rx', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)');
ylabel('Median Error (m)');
title('Median Localization Error');
legend('Golomb Array', 'Evenly Spaced Array');
grid on;

subplot(1, 3, 3);
plot(snr_values, variance_errors_golomb, '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(snr_values, variance_errors_even, '-rx', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)');
ylabel('Error Variance (m^2)');
title('Error Variance');
legend('Golomb Array', 'Evenly Spaced Array');
grid on;

%% Helper Functions
function signals = simulate_propagation(array, source, event_signal, fs, speed_of_sound)
    num_mics = size(array, 1);
    signal_length = length(event_signal);
    max_delay = round(max(vecnorm(array - source, 2, 2)) / speed_of_sound * fs);
    signals = zeros(num_mics, signal_length + max_delay); % Allocate signals with padding
    
    for mic = 1:num_mics
        distance = norm(array(mic, :) - source); % Distance from source to mic
        delay = round(distance / speed_of_sound * fs); % Delay in samples
        signals(mic, delay + (1:signal_length)) = event_signal; % Add delayed signal
    end
    
    signals = signals(:, 1:signal_length); % Trim signals to match length
end

function estimated_position = grid_search_localization(array, signals, radius, event_signal, fs, speed_of_sound)
    x_grid = linspace(-2*radius, 2*radius, 50);
    y_grid = linspace(-2*radius, 2*radius, 50);
    [X, Y] = meshgrid(x_grid, y_grid);
    search_positions = [X(:), Y(:)];
    
    grid_energies = zeros(size(search_positions, 1), 1);
    for i = 1:size(search_positions, 1)
        test_position = search_positions(i, :);
        combined_signal = zeros(size(event_signal));
        
        for mic = 1:size(array, 1)
            distance = norm(array(mic, :) - test_position);
            delay = round(distance / speed_of_sound * fs);
            combined_signal = combined_signal + circshift(signals(mic, :), -delay);
        end
        
        grid_energies(i) = sum(combined_signal.^2);
    end
    
    [~, max_idx] = max(grid_energies);
    estimated_position = search_positions(max_idx, :);
end
%}