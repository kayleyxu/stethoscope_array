clear; close all;

%% Set Random Seed for Reproducibility
seed = 7;
rng(seed);
n_sources = 10;
%n_sources = 100;
num_trials = 2; 
snr_values = [-3, 0, 3]; 
%snr_values = [-3, 0, 3, 6];
frequencies = [150, 250, 400];
%frequencies = [150, 250, 400, 600];
results=[];
results_for_ml=[];

%% Define Output Directory
output_dir = 'simulation_results_ml';
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
c = 4000; % 4000cm/s = 40 m/s
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

side_mic = [-12, 2.3, 10];  % Mic on the side (z ≠ 0)

%% Microphone Positions
position_sets{1} = [[costas_array, zeros(size(costas_array, 1), 1)]; side_mic];

errors = zeros(1, n_sources, length(snr_values), length(frequencies));
all_predictions = cell(1, num_trials, n_sources);

%% Simulation Loop
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

    for f_idx = 1:length(frequencies)
        freq = frequencies(f_idx);
        num_zeros = 20; 
        t_sin = t(1:end-num_zeros); 
        source_signal = [zeros(num_zeros,1); 100 * sin(2 * pi * freq * t_sin)];            
        
        for snr_idx = 1:length(snr_values)
            snr = snr_values(snr_idx);
            trial_errors = zeros(num_trials, 1);
            
            for t_idx = 1:num_trials
                n_mics = size(position_sets{1}, 1);
                mic_signals = zeros(length(t), n_mics);
                delays = zeros(n_mics, 1);
                
                for mic = 1:n_mics
                    distance = norm(position_sets{1}(mic, :) - source_position);
                    delay = distance / c; 
                    delays(mic) = delay;
                    attenuation = 1;
                    mic_signals(:, mic) = attenuation * circshift(source_signal, round(delay * fs));
                end

                mic_signals = awgn(mic_signals, snr, 'measured');
                %% Flattened Raw Microphone Signals
                raw_mic_signals = reshape(mic_signals', 1, []); 
                %% SRP-DAS Localization
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
                        distance = norm(position_sets{1}(mic, :) - test_position);
                        delay = round(distance / c * fs);
                        combined_signal = combined_signal + circshift(mic_signals(:, mic), -delay);
                    end
                    grid_energies(i) = sum(combined_signal.^2);
                end

                [~, max_idx] = max(grid_energies);
                predicted_position_srp = search_positions(max_idx, :);

                %% Step 2: Estimate Delays Using Cross-Correlation
                estimated_delays = zeros(n_mics, 1);
                ref_signal = mic_signals(:, 1);
                for mic = 1:n_mics
                    [xc, lags] = xcorr(ref_signal, mic_signals(:, mic));
                    [~, max_idx] = max(xc);
                    estimated_delays(mic) = -lags(max_idx) / fs;
                end
           
                %% Step 3: Align Delays to Closest Mic
                [~, closest_mic_idx] = min(estimated_delays);
                expected_delay = norm(position_sets{1}(closest_mic_idx, :) - predicted_position_srp) / c;
                delay_offset = expected_delay - estimated_delays(closest_mic_idx);
                adjusted_delays = estimated_delays + delay_offset;
                
%{
delays = adjusted_delays;

                %% Triangulation-Based Localization
                [~, sorted_idx] = sort(delays(1:5));  
                selected_mics_idx = [sorted_idx(1:3);n_mics];  
                selected_mics = position_sets{1}(selected_mics_idx, :);

                ref_mic = selected_mics(1, :);
                A = zeros(3, 3);
                b = zeros(3, 1);

                for i = 2:4
                    mic_i = selected_mics(i, :);
                    tdoa = delays(selected_mics_idx(i)) - delays(selected_mics_idx(1));
                    d_ref = c * delays(selected_mics_idx(1)); 
                    d_i = c * delays(selected_mics_idx(i)); 
                    A(i-1, :) = 2 * (mic_i - ref_mic);
                    b(i-1) = (d_ref^2 - d_i^2) - norm(ref_mic)^2 + norm(mic_i)^2;
                end

                predicted_position_tri = A \ b;
                predicted_position_tri = predicted_position_tri(:)';

                %% Averaged Position
                predicted_position_combined = (predicted_position_srp + predicted_position_tri) / 2;

                %% Error Computation
                error_srp = norm(predicted_position_srp - source_position);
                error_tri = norm(predicted_position_tri - source_position);
                error_combined = norm(predicted_position_combined - source_position);

                fprintf('Source: (%.4f, %.4f, %.4f)\n', source_position);
                fprintf('SRP-DAS: (%.4f, %.4f, %.4f), Error: %.4f\n', predicted_position_srp, error_srp);
                fprintf('Triangulation: (%.4f, %.4f, %.4f), Error: %.4f, rcond(A): %.4f \n', predicted_position_tri, error_tri, rcond(A));
                fprintf('Combined: (%.4f, %.4f, %.4f), Error: %.4f\n\n', predicted_position_combined, error_combined);
 distance_srp_tri = norm(predicted_position_srp - predicted_position_tri); % Compute distance between SRP and Triangulation

%results = [results; source_position, predicted_position_srp, error_srp, ...
%           predicted_position_tri(:)', error_tri, rcond(A), predicted_position_combined, error_combined, distance_srp_tri];
%}
                %trial_errors(t_idx) = error_combined;
                %all_predictions{1, t_idx, s} = predicted_position_combined;
                results_for_ml = [results_for_ml; source_position, snr, freq, c, predicted_position_srp, adjusted_delays', reshape(position_sets{1}, 1, []),raw_mic_signals];

            end
            errors(1, s, snr_idx, f_idx) = mean(trial_errors);
        end
    end
end
save(fullfile(output_dir, 'results_for_ml.mat'), 'results_for_ml');

% costas_mic_array_src_trial_display(errors, frequencies, snr_values, {'Costas Array'}, position_sets, source_positions, all_predictions, output_dir)

diary off;
