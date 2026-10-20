%function [cropped_audio_data_20ms, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path, start_time, end_time)
% find out sound speed through media with two mics
% sound speed Processes microphone audio files to find and crop the most active 20ms window.
% 
% Inputs:
%   folder_path  - Path to the folder containing the 16 mic audio files in .raw format.
%   start_time   - (Optional) Start time of the first-order cropping in seconds.
%   end_time     - (Optional) End time of the first-order cropping in seconds.
% 
% Outputs:
%   cropped_audio_data_20ms - A cell array containing the cropped 20ms audio data for all microphones.
%   common_window_time_ms   - The start time of the 20ms window in milliseconds, relative to the original file.
%   latencies_ms            - Latencies for each mic relative to the earliest mic in ms.
%
% Figures:
%   - Waveforms and delays for all mics (50ms window).
%   - STFT spectrogram for each microphone (50ms window).

%% Constants
clear;
%folder_path = '../data/sound_speed_2/'
folder_path = '../data/sound_speed_3/'

num_mics = 2;
%Fs = 64000;            % Sampling rate in Hz
Fs=44100
bits_per_sample = 16;  % 16-bit audio
sample_format = 'int16'; % Format of raw data (signed 16-bit integers)
% these calibration values are obtained by running source at center scenario
cal_samples=[0     0];
%cal_samples=[0     0     0     0     0     0     1     1     2     1     1     1     1     1     0     0];
%cal_ms=[0         0         0         0         0         0    0.0208    0.0208    0.0417    0.0208    0.0208    0.0208    0.0208    0.0208         0         0];

%% Step 1: Load and first-order crop the 16 audio files
%disp('Loading and optionally cropping audio files...');
audio_files = cell(num_mics, 1);
start_time=5;
end_time=15;

% Crop indices (if start_time and end_time are provided)
if ~isempty(start_time) && ~isempty(end_time)
    crop_start_idx = round(start_time * Fs)
    crop_end_idx = round(end_time * Fs);
else
    crop_start_idx = 1;
    crop_end_idx = Inf;
end

% Load and optionally crop each file
for mic = 1:num_mics
    file_name = sprintf('audio_opt%d.raw', mic-1);
    file_path = fullfile(folder_path, file_name);
    
    % Open and read the .raw file
    fid = fopen(file_path, 'rb');
    if fid == -1
        error(['Unable to open file: ', file_path]);
    end
    raw_data = fread(fid, Inf, sample_format); % Read the entire file
    fclose(fid);
    
    % Normalize audio data to range [-1, 1]
    raw_audio = double(raw_data) / (2^(bits_per_sample - 1));
    
    % Apply first-order cropping
    crop_end_idx = min(crop_end_idx, length(raw_audio)); % Ensure the end index is valid
    audio_files{mic} = raw_audio(crop_start_idx:crop_end_idx);
    %%%%%    cal_ms = round(cal(mic) * Fs/1000); % cal value in ms
    % calibration compensation
    audio_files{mic} = circshift(audio_files{mic},-cal_samples(mic));
end
    original_length = length(audio_files{1});
%{
    if nargin > 1
    disp(['First-order cropping applied: Start Time = ', num2str(start_time), ' s, End Time = ', num2str(end_time), ' s']);
else
    disp('No first-order cropping applied (using full signals).');
end
%}
highpass_cutoff = 50;  % High-pass filter cutoff frequency in Hz (adjust as needed)
filter_order = 4;  % Order of the Butterworth filter

% Design the high-pass filter
[b, a] = butter(filter_order, highpass_cutoff / (Fs / 2), 'high')
% Apply the filter to both mic signals
audio_files{1} = filter(b, a, audio_files{1});
audio_files{2}= filter(b, a, audio_files{2});
    figure
subplot(2,1,1)
plot(audio_files{1})
subplot(2,1,2)
plot(audio_files{2})

%% Step 2: Identify common 20ms window with signal of interst processing
% search with 10ms sliding window
%disp('Identifying the 20ms window of interest with 10ms sliding...');
window_size = round(0.02 * Fs); % 20ms window in samples
sliding_size = round(0.01 * Fs); % 10ms sliding window in samples
global_activity = zeros(1, original_length - sliding_size + 1);

% Compute the sum of squared amplitudes for all microphones
for mic = 2:num_mics
    signal = audio_files{mic};
    local_activity = zeros(1, length(signal) - sliding_size + 1);
    
    % Sum of squared amplitudes within the sliding window
    for i = 1:length(local_activity)
        local_activity(i) = sum(signal(i:i + sliding_size - 1).^2);
    end
    
    % Accumulate global activity across all microphones
    global_activity = global_activity + local_activity;
end

% Find the sliding window with the maximum global activity
[~, max_idx] = max(global_activity);

% Compute common 20ms window indices
activity_peak_global = max_idx + round(sliding_size / 2);
start_idx_20ms = max(1, activity_peak_global - floor(window_size / 2)); % Start of centered window
end_idx_20ms = min(original_length, start_idx_20ms + window_size - 1); % End of 20ms window

% Extend the window to 50ms for plotting and STFT
extension_50ms = round(0.025 * Fs); % 25ms extension on each side
start_idx_50ms = max(1, start_idx_20ms - extension_50ms);
end_idx_50ms = min(original_length, end_idx_20ms + extension_50ms);

cropped_audio_files_50ms = cell(num_mics, 1);
cropped_audio_data_20ms = cell(num_mics, 1);
for mic = 1:num_mics
    cropped_audio_files_50ms{mic} = audio_files{mic}(start_idx_50ms:end_idx_50ms);
    cropped_audio_data_20ms{mic} = audio_files{mic}(start_idx_20ms:end_idx_20ms);
end

% Display the time of the common 20ms window
common_window_time_ms = ((start_idx_20ms - 1) + crop_start_idx - 1) / Fs * 1000; % Convert to ms
%disp(['Common 20ms window identified: Start Time = ', num2str(round(common_window_time_ms, 2)), ' ms']);

%% Step 3: Compute latencies using cross-correlation (20ms window)
%disp('Computing latencies using cross-correlation (20ms window)...');
reference_signal = cropped_audio_data_20ms{1};
latencies_ms = zeros(1, num_mics); % Store latencies in ms
correlation_results = cell(1, num_mics); % Store cross-correlation results

for mic = 1:num_mics
    % Extract 20ms window for the current mic
    mic_window=cropped_audio_data_20ms{mic};
    % Compute cross-correlation using the 20ms windows
    [correlation, lag] = xcorr(mic_window, reference_signal, 'coeff'); % Normalized correlation
    correlation_results{mic} = struct('correlation', correlation, 'lag', lag); % Store results
    
    % Find the lag with the maximum correlation
    [~, max_idx] = max(correlation);
    latencies_samples(mic) = lag(max_idx);
    latencies_ms(mic) = latencies_samples(mic) / Fs * 1000; % Convert lag to ms
 end
% Adjust latencies relative to the earliest arriving mic
[~, earliest_mic_idx] = min(latencies_samples);
latencies_samples = latencies_samples - latencies_samples(earliest_mic_idx);
latencies_ms = latencies_ms - latencies_ms(earliest_mic_idx);
%d=0.0638
%d=0.073
d=0.078

max(latencies_ms)
c=d/(max(latencies_ms)/1000)
% Display correlation matrix results for debugging
%{
disp('Cross-correlation results (debug):');
for mic = 1:num_mics
    fprintf('Mic %d: Max correlation = %.4f, Lag = %.4f ms\n', ...
            mic-1, max(correlation_results{mic}.correlation), ...
            correlation_results{mic}.lag(find(correlation_results{mic}.correlation == max(correlation_results{mic}.correlation), 1)) / Fs * 1000);
    latencies_ms(mic)
    latencies_samples(mic)
end
%}
%% Step 4: Plot waveforms and delays
%disp('Plotting waveforms and delays...');
figure;

hold on;
time_vector = (0:length(cropped_audio_files_50ms{1}) - 1) / Fs * 1000; % Time vector for the 50ms window in ms

% Calculate maximum absolute amplitude across all waveforms for scaling
max_amplitude = max(cellfun(@(x) max(abs(x)), cropped_audio_files_50ms));

for mic = 1:num_mics
    % Offset for clarity
    offset = (mic - 1) * (2 * max_amplitude); % Scale offset based on maximum amplitude
    
    % Plot waveform
    plot(time_vector, cropped_audio_files_50ms{mic} + offset, 'DisplayName', ['Mic ', num2str(mic-1)]);
    
    % Annotate latency
    if mic == earliest_mic_idx
        %annotation_text = sprintf('%d samples, %.4f ms (Earliest)', latencies_samples(mic), latencies_ms(mic));
        annotation_text = sprintf('%.4f ms (Earliest)', latencies_ms(mic));
    else
        %annotation_text = sprintf('%d samples, %.4f ms', latencies_samples(mic),latencies_ms(mic));
        annotation_text = sprintf('%.4f ms',latencies_ms(mic));
    end
    text(time_vector(1), offset, annotation_text, 'Color', 'red', 'FontSize', 10, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        % Annotate mic index at the right of the plot
    text(time_vector(end), offset, sprintf('Mic %d', mic-1), 'Color', 'blue', 'FontSize', 10, ...
         'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
end

% Finalize plot
title('Waveform and Arrival Latency');
xlabel('Time (ms)');
ylabel('Amplitude');


%% Step 5: Frequency Analysis Using STFT
%disp('Analyzing the frequency of the crackle sound using STFT...');
fft_size = 64; % STFT FFT size for ~1ms resolution
window_size = fft_size; % Match the FFT size
overlap_size = round(fft_size / 2); % 50% overlap

% Plot STFT spectrogram for all microphones
figure;
sgtitle('STFT Spectrogram for Each Microphone (50ms Window)');
for mic = 1:num_mics
    % Perform STFT on the 50ms window signal
    signal = cropped_audio_files_50ms{mic};
    [S, F, T] = stft(signal, Fs, 'Window', hamming(window_size), ...
                     'OverlapLength', overlap_size, 'FFTLength', fft_size);
    
    % Magnitude of the STFT
    magnitude = abs(S);
    
    % Plot the spectrogram
    subplot(4, 4, mic);
    surf(T * 1000, F, 20 * log10(magnitude), 'EdgeColor', 'none'); % Plot in ms
    axis tight;
    view(0, 90);
    colormap jet;
    colorbar;
    title(['Mic ', num2str(mic-1)]);
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    ylim([0, Fs/2]); % Show only positive frequencies
end
figure
subplot(2,1,1)
plot(audio_files{1})
subplot(2,1,2)
plot(audio_files{2})