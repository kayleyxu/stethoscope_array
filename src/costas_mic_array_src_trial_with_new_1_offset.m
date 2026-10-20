clear; close all;

%% Set Random Seed for Reproducibility
seed = 7;
rng(seed);
n_sources = 10;
num_trials = 2; 
snr_values = [-3, 0, 3];
frequencies = [150, 250, 400];
num_additional_arrangements = 10; 

%% Define Output Directory
output_dir = 'simulation_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
diary(fullfile(output_dir, 'simulation_log.txt'));
diary on;

%% Acoustic Box Setup
box_length = 25; 
box_width = 35; 
box_height = 20; 

%% Simulation Parameters
c = 4000; %cm/s
fs = 44100; 
t = (0:1/fs:0.01)'; 

%% Microphone Placement Configurations
scaling_factor = 3;
x_offset = -10;
y_offset = -4;

costas_array = [
     1, 1;
     2, 4;
     3, 2;
     4, 3;
     5, 5
] * scaling_factor + [x_offset, y_offset];

%% Original 10 Perturbed Arrays (±0.5)
perturbations = round(1 * rand(num_additional_arrangements, size(costas_array, 1), 2)) - 0.5;
additional_arrays = cell(num_additional_arrangements, 1);

for i = 1:num_additional_arrangements
    additional_arrays{i} = costas_array + squeeze(perturbations(i, :, :));
end

%% Placements and Position Sets (Original 11 Arrays)
placements = [{'Costas Array A0'}, arrayfun(@(i) sprintf('Perturbed A%d', i), 1:num_additional_arrangements, 'UniformOutput', false)];
position_sets = cell(1, num_additional_arrangements + 1);
% 5 mics + 6th on side
position_sets{1} = [[costas_array, zeros(size(costas_array, 1), 1)]; [-12, 2.3, 10]];
% 5 mics
% position_sets{1} = [costas_array, zeros(size(costas_array, 1), 1)]; 
for i = 1:num_additional_arrangements
    % 5 mics + 6th on side
    position_sets{i + 1} = [[additional_arrays{i}, zeros(size(additional_arrays{i}, 1), 1)]; [-12, 2.3, 10]];
    % 5 mics
    % position_sets{i + 1} = [additional_arrays{i}, zeros(size(additional_arrays{i}, 1), 1)];
end

errors = zeros(length(position_sets), n_sources, length(snr_values), length(frequencies));
all_predictions = cell(length(position_sets), num_trials, n_sources);

%% Save RNG State Before Original Trials
rng_state_before_trials = rng;

%% Simulation Loop for Original Arrays
for s = 1:n_sources
    z = rand(1) * 8 + 4;
    k = convhull(costas_array(:, 1), costas_array(:, 2));
    
    while true
        x = min(costas_array(:, 1)) + rand(1) * (max(costas_array(:, 1)) - min(costas_array(:, 1)));
        y = min(costas_array(:, 2)) + rand(1) * (max(costas_array(:, 2)) - min(costas_array(:, 2)));
        if inpolygon(x, y, costas_array(k, 1), costas_array(k, 2))
            break;
        end
    end
    
    source_position = [x, y, z];
    source_positions(s,:) = source_position;

    for p = 1:length(position_sets)
        array_positions = position_sets{p};
        
        for f_idx = 1:length(frequencies)
            freq = frequencies(f_idx);
            source_signal = 100 * sin(2 * pi * freq * t);
            for snr_idx = 1:length(snr_values)
                snr = snr_values(snr_idx);
                trial_errors = zeros(num_trials, 1);
                
                for t_idx = 1:num_trials
                    n_mics = size(array_positions, 1);
                    mic_signals = zeros(length(t), n_mics);
                    
                    for mic = 1:n_mics
                        distance = norm(array_positions(mic, :) - source_position);
                        delay = round(distance / c * fs);
                        attenuation = 1; 
                        mic_signals(:, mic) = attenuation * circshift(source_signal, delay);
                    end

                    mic_signals = awgn(mic_signals, snr, 'measured');

                    x_grid = linspace(-10, 10, 40);
                    y_grid = linspace(-10, 10, 40);
                    z_grid = linspace(5, 12, 14);

                    [X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);
                    search_positions = [X(:), Y(:), Z(:)];

                    grid_energies = zeros(size(search_positions, 1), 1);
                    for i = 1:size(search_positions, 1)
                        test_position = search_positions(i, :);
                        combined_signal = zeros(size(t));

                        for mic = 1:n_mics
                            distance = norm(array_positions(mic, :) - test_position);
                            delay = round(distance / c * fs);
                            combined_signal = combined_signal + circshift(mic_signals(:, mic), -delay);
                        end
                        grid_energies(i) = sum(combined_signal.^2);
                    end

                    [~, max_idx] = max(grid_energies);
                    predicted_position = search_positions(max_idx, :);

                    trial_errors(t_idx) = norm(predicted_position - source_position);
                    all_predictions{p, t_idx, s} = predicted_position;
                end
                errors(p, s, snr_idx, f_idx) = mean(trial_errors);
            end
        end
    end
end

%% Add 10 New Additional Arrays (±1 Perturbations)
rng_state = rng; % Save RNG state before generating new perturbations
new_perturbations = 2* round(rand(10, size(costas_array, 1), 2)) - 1;

for i = 1:10
    additional_arrays{num_additional_arrangements + i} = costas_array + squeeze(new_perturbations(i, :, :));
end
rng(rng_state_before_trials); % Restore RNG state to ensure original trials are unchanged

%% Append New Arrays
for i = 1:10
    %placements{end+1} = sprintf('New Perturbed Array %d', i);
    placements{end+1} = sprintf('Perturbed A%d', i+10);
    % 5 mics + 6th on the side
    position_sets{end+1} = [[additional_arrays{num_additional_arrangements + i}, zeros(size(costas_array, 1), 1)]; [-12, 2.3, 10]];
    % 5 mics
    %position_sets{end+1} = [additional_arrays{num_additional_arrangements + i}, zeros(size(costas_array, 1), 1)];
    end

%% Simulation for New Arrays Only
for s = 1:n_sources
    source_position = source_positions(s, :);

    for p = (num_additional_arrangements + 2):(num_additional_arrangements + 11)
        array_positions = position_sets{p};

        for f_idx = 1:length(frequencies)
            freq = frequencies(f_idx);
            source_signal = 100 * sin(2 * pi * freq * t);

            for snr_idx = 1:length(snr_values)
                snr = snr_values(snr_idx);
                trial_errors = zeros(num_trials, 1);

                for t_idx = 1:num_trials
                    n_mics = size(array_positions, 1);
                    mic_signals = zeros(length(t), n_mics);

                    for mic = 1:n_mics
                        distance = norm(array_positions(mic, :) - source_position);
                        delay = round(distance / c * fs);
                        attenuation = 1;
                        mic_signals(:, mic) = attenuation * circshift(source_signal, delay);
                    end

                    mic_signals = awgn(mic_signals, snr, 'measured');

                    x_grid = linspace(-10, 10, 40);
                    y_grid = linspace(-10, 10, 40);
                    z_grid = linspace(5, 12, 14);

                    [X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);
                    search_positions = [X(:), Y(:), Z(:)];

                    grid_energies = zeros(size(search_positions, 1), 1);
                    for i = 1:size(search_positions, 1)
                        test_position = search_positions(i, :);
                        combined_signal = zeros(size(t));

                        for mic = 1:n_mics
                            distance = norm(array_positions(mic, :) - test_position);
                            delay = round(distance / c * fs);
                            combined_signal = combined_signal + circshift(mic_signals(:, mic), -delay);
                        end
                        grid_energies(i) = sum(combined_signal.^2);
                    end

                    [~, max_idx] = max(grid_energies);
                    predicted_position = search_positions(max_idx, :);

                    trial_errors(t_idx) = norm(predicted_position - source_position);
                    all_predictions{p, t_idx, s} = predicted_position;
                end
                errors(p, s, snr_idx, f_idx) = mean(trial_errors);
            end
        end
    end
end

costas_mic_array_src_trial_display(errors, frequencies, snr_values, placements, position_sets, source_positions, all_predictions, output_dir)

diary off;
