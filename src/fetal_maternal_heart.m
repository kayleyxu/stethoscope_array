% Use differential beamforming for spatial nulling
clear; clc;

% Parameters
Fs = 44100; % Sampling rate in Hz
c = 40;            % Speed of sound (m/s)
duration = 10;
volume_scale = 10; % Volume scaling factor         
duration = 10; % Total duration in seconds (2 breathing cycles)

% Array Parameters
N = 2;                  % Number of microphones
ang_n = -90; %89.9 % Null angle in degrees
%dd = 0.016;%0.032; % physical spacing between two mics
dd=0.005;
micarrayd = phased.ULA(N, dd);

% Frequency range of interest
f_min = 1;
f_max = 1000;
frequencies = f_min:f_max;

% source location
%source_1 = [0, -0.04, 0.05]; % Source at -90° (heart)
%source_2 = [0.03, 0.06, 0.05];  % Source at +90° (lung)
%source_1 = [0, -0.07, 0.05]; % Source at -90° (heart)
%source_2 = [0.03, 0.06, 0.05];  % Source at +90° (lung)
%source_1 = [0, -0.17, 0.05]; % Source at -90° (mother heart)
%source_2 = [0.03, 0.13, 0.05];  % Source at +90° (fetus heart)
source_1 = [0, -0.20, 0.04]; % Source at -90° (mother heart)
source_2 = [0.03, 0.05, 0.04];  % Source at +90° (fetus heart)
%{
% crackle sound recorded from in-vitro experiment
folder_path = '..\data\2024_11_26_23_51_25_single_crackle';
[cropped_audio_data, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path);
% take from the first mic
crackle_sound = cropped_audio_data{1};

t = 0:1/Fs:duration; % Time vector
disp('Generating synthetic breathing and heartbeat sounds...');

% Breathing cycle parameters
breath_rate = 2 / duration; % 2 breaths in total duration
cycle_duration = 1 / breath_rate; % Total duration of one cycle
inhale_ratio = 0.4; % Inhale duration ratio (40% inhale, 60% exhale)
inhale_duration = inhale_ratio * cycle_duration; % Inhale duration
exhale_duration = (1 - inhale_ratio) * cycle_duration; % Exhale duration

% Low-pass filter parameters for realistic breathing sound
fc_low = 500; % Cutoff frequency (Hz) for breathing sound
[b_low, a_low] = butter(4, fc_low / (Fs / 2), 'low'); % Low-pass filter for breathing

% Parameters for inhale and exhale
inhale_amplitude = 0.4; % Softer inhale
exhale_amplitude = 0.8; % Stronger exhale

% Generate breathing sound with crackles
breathing_with_crackles = zeros(size(t));
crackle_times = []; % Store times where crackles are added
for cycle = 1:2 % Two breathing cycles
    % Inhale phase
    cycle_start_time = (cycle - 1) * cycle_duration;
    inhale_end_time = cycle_start_time + inhale_duration;
    inhale_indices = (t >= cycle_start_time) & (t < inhale_end_time);
    inhale_samples = sum(inhale_indices);
    inhale_time = linspace(0, inhale_duration, inhale_samples);

    % Inhale: Turbulent noise with smooth amplitude modulation
    inhale_wave = inhale_amplitude * (1 - cos(2 * pi * inhale_time / inhale_duration)) .* ...
                  randn(1, inhale_samples); % White noise modulated for inhale
    inhale_filtered = filter(b_low, a_low, inhale_wave); % Apply low-pass filter for stethoscope effect
    breathing_with_crackles(inhale_indices) = inhale_filtered;

    % Exhale phase
    exhale_start_time = inhale_end_time;
    exhale_end_time = cycle_start_time + cycle_duration;
    exhale_indices = (t >= exhale_start_time) & (t < exhale_end_time);
    exhale_samples = sum(exhale_indices);
    exhale_time = linspace(0, exhale_duration, exhale_samples);

    % Exhale: Harsher turbulent noise with decaying amplitude
    exhale_wave = exhale_amplitude * exp(-exhale_time / (exhale_duration / 3)) .* ...
                  randn(1, exhale_samples); % White noise modulated for exhale

    % Apply a Hanning window to the first 10% of the exhale to smooth the onset
    hanning_length = round(0.1 * exhale_samples); % First 10% of exhale samples
    hanning_window = hann(2 * hanning_length)'; % Create a Hanning window
    hanning_window = hanning_window(1:hanning_length); % Use the first half of the window

    % Apply the Hanning window to the start of the exhale_wave
    exhale_wave(1:hanning_length) = exhale_wave(1:hanning_length) .* hanning_window;

    exhale_filtered = filter(b_low, a_low, exhale_wave); % Apply low-pass filter for stethoscope effect
    breathing_with_crackles(exhale_indices) = exhale_filtered;

    % Add real crackle sound to simulated breathing sound
    crackle_start_time = exhale_start_time + (exhale_duration / 2); % Mid-exhale
    crackle_start_idx = round(crackle_start_time * Fs);
    crackle_end_idx = crackle_start_idx + length(crackle_sound) - 1;
    if crackle_end_idx <= length(t)
        aligned_crackle = 15 * crackle_sound(1:(crackle_end_idx - crackle_start_idx + 1));
        breathing_with_crackles(crackle_start_idx:crackle_end_idx) = ...
            breathing_with_crackles(crackle_start_idx:crackle_end_idx) + aligned_crackle';
        crackle_times = [crackle_times, crackle_start_time]; % Store time in seconds
    end
end

% Normalize the sound
breathing_with_crackles = breathing_with_crackles / max(abs(breathing_with_crackles));

% Heartbeat Signal Parameters
heartbeat_rate = 80 / 60;   % Heartbeat frequency in Hz
s1_duration = 0.1; %0.1;
s2_duration = 0.07;  
systolic_interval = 0.3; %0.3;
diastolic_interval = (1 / heartbeat_rate) - (s1_duration + s2_duration + systolic_interval);

% Generate Heartbeat Signals
heartbeat_pattern = [ones(1, round(s1_duration * Fs)), ...
                         zeros(1, round(systolic_interval * Fs)), ...
                         ones(1, round(s2_duration * Fs)), ...
                         zeros(1, round(diastolic_interval * Fs))];
heartbeat_signal = repmat(heartbeat_pattern, 1, ceil(length(t) / length(heartbeat_pattern)));
heartbeat_signal = circshift(heartbeat_signal(1:length(t)), round(0.1 * Fs));
%heartbeat_signal = bandpass(heartbeat_signal, [50, 80], Fs); 
heartbeat_signal = bandpass(heartbeat_signal, [50, 200], Fs);
heartbeat_signal = 0.5*heartbeat_signal/ max(abs(heartbeat_signal));

%}

% Heartbeat Signal Parameters
t = 0:1/Fs:duration; % Time vector
mom_heartbeat_rate = 80 / 60;   % Heartbeat frequency in Hz for mom
fetus_heartbeat_rate = 140 / 60; % Heartbeat frequency in Hz for fetus
s1_duration_mom = 0.1;
s2_duration_mom = 0.07;
s1_duration_fetus = 0.06; % Duration of S1 sound (seconds)
s2_duration_fetus = 0.04; % Duration of S2 sound (seconds)     
systolic_interval_mom = 0.3;
diastolic_interval_mom = (1 / mom_heartbeat_rate) - (s1_duration_mom + s2_duration_mom + systolic_interval_mom);
systolic_interval_fetus = 0.2;
diastolic_interval_fetus = (1 / fetus_heartbeat_rate) - (s1_duration_fetus + s2_duration_fetus + systolic_interval_fetus);   % Silence between beats (seconds)

% Generate Heartbeat Signals
heartbeat_pattern_mom = [ones(1, round(s1_duration_mom * Fs)), ...
                         zeros(1, round(systolic_interval_mom * Fs)), ...
                         ones(1, round(s2_duration_mom * Fs)), ...
                         zeros(1, round(diastolic_interval_mom * Fs))];
heartbeat_signal_mom = repmat(heartbeat_pattern_mom, 1, ceil(length(t) / length(heartbeat_pattern_mom)));
heartbeat_signal_mom = circshift(heartbeat_signal_mom(1:length(t)), round(0.1 * Fs));
heartbeat_signal_mom = bandpass(heartbeat_signal_mom, [50, 200], Fs);

heartbeat_pattern_fetus = [ones(1, round(s1_duration_fetus * Fs)), ...
                           zeros(1, round(systolic_interval_fetus * Fs)), ...
                           ones(1, round(s2_duration_fetus * Fs)), ...
                           zeros(1, round(diastolic_interval_fetus * Fs))];
heartbeat_signal_fetus = repmat(heartbeat_pattern_fetus, 1, ceil(length(t) / length(heartbeat_pattern_fetus)));
heartbeat_signal_fetus = circshift(heartbeat_signal_fetus(1:length(t)), round(0.2 * Fs));
heartbeat_signal_fetus = bandpass(heartbeat_signal_fetus, [50, 200], Fs);

breathing_with_crackles = heartbeat_signal_fetus/ max(abs(heartbeat_signal_fetus));
heartbeat_signal=heartbeat_signal_mom/ max(abs(heartbeat_signal_mom));

% just dummy lambda to get pos', ignore dummy_w
lambda = c / 70; 
dd_f = dd / lambda;
[dummy_w, pos] = diffbfweights(N, dd_f, ang_n, 'ArrayGeometry', 'ULA'); % Differential weight
mic_pos = pos' * lambda; % mic_pos is Nx3

delays_1 = vecnorm(mic_pos - source_1, 2, 2) / c; % Delays for heart (-90°)
delays_2 = vecnorm(mic_pos - source_2, 2, 2) / c; % Delays for lung (+90°)

% Simulate Microphone Signals at the two mics
mic_signals = zeros(N, length(t));
for i = 1:N
    mic_signals(i, :) = circshift(heartbeat_signal, round(delays_1(i) * Fs)) + ...
                        circshift(breathing_with_crackles, round(delays_2(i) * Fs));
end

% Apply the Hilbert Transform to Extend Signals to Complex Domain
%analytic_signals = hilbert(mic_signals')'; % NxT

% STFT parameters
%window_length = 1024;
%overlap = 512;
%nfft = 1024;
window_length = 1024*4;
overlap = 512*4;
nfft = 1024*4;
% Perform STFT on microphone signals
[S1, f, time] = stft(mic_signals(1, :), Fs, ...
    'Window', hamming(window_length, 'periodic'), ...
    'OverlapLength', overlap, ...
    'FFTLength', nfft, ...
    'Centered', false);  % Ensure zero frequency is at the start
[S2, ~, ~] = stft(mic_signals(2, :), Fs, ...
    'Window', hamming(window_length, 'periodic'), ...
    'OverlapLength', overlap, ...
    'FFTLength', nfft, ...
    'Centered', false);  % Ensure zero frequency is at the start

% Initialize beamformed STFT
S_beamformed = zeros(size(S1));
beamformed_output = zeros(size(mic_signals(1,:)));

% Ensure f_max is within Nyquist frequency
f_max = min(f_max, Fs / 2);

% Compute weights for each frequency bin in STFT
weights = zeros(length(f),N);
for k = 1:length(f)
    if f(k) >= f_min && f(k) <= f_max
        lambda = c / f(k);             % Wavelength for current frequency
        dd_f = dd / lambda;            % Adjusted spacing
        W = diffbfweights(N, dd_f, ang_n, 'ArrayGeometry', 'ULA');  % Dynamic weights
        weights(k,:)=W';
        S_beamformed(k, :) = S1(k, :) * W(1)'+ S2(k, :) * W(2)';
    else
        % Apply attenuation for frequencies outside the target range
        weights(k,:)=[0.1;0.1]';
        S_beamformed(k, :) = 0.1 * (S1(k, :) + S2(k, :));
    end
end

% Perform inverse STFT to obtain time-domain beamformed signal
%beamformed_output = istft(S_beamformed, Fs, 'Window', hamming(window_length, 'periodic'), 'OverlapLength', overlap, 'FFTLength', nfft);
beamformed_output = istft(S_beamformed, Fs, ...
    'Window', hamming(window_length, 'periodic'), ...
    'OverlapLength', overlap, ...
    'FFTLength', nfft, ...
    'Centered', false);

% Normalize the beamformed output
beamformed_output = real(beamformed_output);
beamformed_output = beamformed_output / max(abs(beamformed_output));

% Play the beamformed output
%sound(beamformed_output, Fs);

% Plot Mic and Source Locations in 3D
figure;
subplot(1, 2, 1);
hold on;
plot3(source_2(1), source_2(2), source_2(3), 'cd',  'MarkerSize', 16,'MarkerFaceColor','c', 'DisplayName', 'Lung Source');
plot3(source_1(1), source_1(2), source_1(3), 'md', 'MarkerSize', 16, 'MarkerFaceColor','m','DisplayName', 'Heart Source');
plot3(mic_pos(1, 1), mic_pos(1, 2), mic_pos(1, 3), 'ro', 'MarkerSize', 16, 'MarkerFaceColor','r','DisplayName', 'Microphone 1');
plot3(mic_pos(2, 1), mic_pos(2, 2), mic_pos(2, 3), 'bo', 'MarkerSize', 16, 'MarkerFaceColor','b','DisplayName', 'Microphone 2');
hold off;
title('Mic and Source Locations (3D)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
grid on;
legend('show');
view(3);
axis([-.15 .15 -.15 .15 0 0.1]);

% Plot Beam Pattern
subplot(1, 2, 2);
pattern(micarrayd, f(2), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(2,:)');
hold on;
pattern(micarrayd, f(16), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(16,:)');
pattern(micarrayd, f(5*4), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(5*4,:)');
pattern(micarrayd, f(10*4), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(10*4,:)');
pattern(micarrayd, f(15*4), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(15*4,:)');
pattern(micarrayd, f(20*4), -180:180, 0, 'PropagationSpeed', c, 'Weights', weights(20*4,:)');

title('Beam Pattern');

% Plot Signals
figure;
for i = 1:5
    subplot(5, 2, (i-1)*2 + 1);
    switch i
        case 1
            plot(t, breathing_with_crackles,'c');
            title('Lung Sound (Source)');
        case 2
            plot(t, heartbeat_signal,'m');
            title('Heartbeat Signal (Source)');
            ylim([-1 1]);
        case 3
            plot(t, mic_signals(1, :),'r');
            title('Microphone 1 (Before Beamforming)');
        case 4
            plot(t, mic_signals(2, :),'b');
            title('Microphone 2 (Before Beamforming)');
        case 5
            %plot(t, beamformed_output);
            plot((0:length(beamformed_output)-1)/Fs, beamformed_output,'g');
            title('Beamformed Output (Post Beamforming)');
    end
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    subplot(5, 2, (i-1)*2 + 2);
    % Enhanced Spectrogram
    window_length = 1024;  % Increased window size for better resolution
    overlap = 256;         % Moderate overlap to highlight transient events
    nfft = 1024;           % Higher FFT size for better frequency resolution

    switch i
        case 1
            spectrogram(breathing_with_crackles, window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Lung Source');
        case 2
            spectrogram(heartbeat_signal, window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Heartbeat Source');
        case 3
            spectrogram(mic_signals(1, :), window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Microphone 1');
        case 4
            spectrogram(mic_signals(2, :), window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Microphone 2');
        case 5
            spectrogram(beamformed_output, window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Beamformed Output');
    end

    % Adjust color limits for better contrast
    caxis([-120 -20]); % Adjust based on observed spectrogram

    % Apply a different colormap
    colormap(jet); 

    ylim([0 3]); % Adjust frequency range if needed
end
%}

% Play Sounds
play_start= 2;
play_end =5.5;
play_duration=play_end - play_start;
disp('Maternal heartbeat sound alone...');
sound(heartbeat_signal(play_start*Fs:play_end*Fs) * volume_scale, Fs);
audiowrite('../data/outputs/sound/demo2_maternal_heart.wav',heartbeat_signal(play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

disp('Fetal heartbeat sound alone...');
sound(breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
audiowrite('../data/outputs/sound/demo2_fetal_heart.wav',breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

disp('Sound at Mic 1 (before beamforming)...');
sound(mic_signals(1, play_start*Fs:play_end*Fs) * volume_scale, Fs);
audiowrite('../data/outputs/sound/demo2_mic1_input.wav', mic_signals(1, play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

%disp('Sound at Mic 2 (before beamforming)...');
%sound(mic_signals(2, play_start*Fs:play_end*Fs)) * volume_scale, Fs);
%audiowrite('../data/outputs/sound/mic2_input.wav', mic_signals(2, play_start*Fs:play_end*Fs) * volume_scale, Fs);
%pause(play_duration + 1);

disp('Beamformed output (after beamforming)...');
sound(beamformed_output(play_start*Fs:play_end*Fs)* volume_scale, Fs);
audiowrite('../data/outputs/sound/demo2_beamformed_output.wav', beamformed_output(play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

disp('Compare with fetus heartbeat sound alone...');
sound(breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);
%{
% Define the target duration for the .wav file
target_duration = 60; % Target duration in seconds

% Extend mom's heartbeat signal to match the target duration
heartbeat_signal_extended = repmat(heartbeat_signal, 1, ceil(target_duration / duration));
heartbeat_signal_extended = heartbeat_signal_extended(1:target_duration * Fs);

% Extend fetus's heartbeat signal to match the target duration
breathing_with_crackles_extended = repmat(breathing_with_crackles, 1, ceil(target_duration / duration));
breathing_with_crackles_extended = breathing_with_crackles_extended(1:target_duration * Fs);

% Save the signals as .wav files
audiowrite('heartbeat_signal.wav', heartbeat_signal_extended * volume_scale, Fs);
audiowrite('breathing_with_crackles.wav', breathing_with_crackles_extended * volume_scale, Fs);

disp('Heartbeat signals saved as .wav files.');
target_duration = 3;
mic1_signals_extended = repmat(mic_signals(1, :), 1, ceil(target_duration / duration));
mic1_signals_extended = mic1_signals_extended(1:target_duration * Fs);
mic2_signals_extended = repmat(mic_signals(2, :), 1, ceil(target_duration / duration));
mic2_signals_extended = mic2_signals_extended(1:target_duration * Fs);
audiowrite('mic1_signals.wav', mic1_signals_extended, Fs);
audiowrite('mic2_signals.wav', mic2_signals_extended, Fs);

figure;
subplot(2,2,1);
pattern(micarrayd, f, -180:180, 0, 'PropagationSpeed', c, 'Weights', w);
subplot(2,2,2);
pattern(micarrayd, f, -180:180, 0, 'PropagationSpeed', c+10, 'Weights', w);
subplot(2,2,3);
pattern(micarrayd, f+10, -180:180, 0, 'PropagationSpeed', c, 'Weights', w);
subplot(2,2,4);
pattern(micarrayd, f+10, -180:180, 0, 'PropagationSpeed', c+20, 'Weights', w);

%}

