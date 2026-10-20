% Use differential beamforming for spatial nulling
clear; clc;

% Parameters
Fs = 44100; % Sampling rate in Hz
volume_scale = 1; % Volume scaling factor         
duration = 4;%5%10; % Total duration in seconds (2 breathing cycles)

% Array Parameters
N = 2;                  % Number of microphones
% Frequency range of interest
f_min = 1;
f_max = 3000;
frequencies = f_min:f_max;
t = 0:1/Fs:duration; % Time vector
mic_signals = zeros(N, length(t));

% in-vitro test code section
%mic1_file = '../data/diff_invitro_2025-02-17_14_13_56/audio_opt1.raw';
%mic2_file = '../data/diff_invitro_2025-02-17_14_13_56/audio_opt0.raw';
% invitro_intf_supress_2: good volume, 2 revb 0.012, mic1= audio_opt8,
% mic2=audio_opt7, start_time =7

% invitro_intf_supress_3: smaller heart volume, 2 revb 0.012
% invitro_intf_supress_4: smaller heart volume, 1 reva+1revB MP 0.005, mic1= audio_opt7,
% mic2=audio_opt10, start_time=14
% invitro_intf_supress_5: 1 reva+1revB CMM 0.005, mic1= audio_opt7,
% mic2=audio_opt11, start_time=? 
% invitro_intf_supress_6: large  heart volume, 1 reva+1revBCMM 0.005, mic1= audio_opt7,
% mic2=audio_opt11, start_time=? duration: 10



%mic1_file = '../data/invitro_intf_supress/audio_opt7.raw';
%mic2_file = '../data/invitro_intf_supress/audio_opt8.raw';
mic1_file = '../data/invitro_intf_supress_2/audio_opt8.raw';
mic2_file = '../data/invitro_intf_supress_2/audio_opt7.raw';
%mic1_file = '../data/invitro_intf_supress_4/audio_opt7.raw';
%mic2_file = '../data/invitro_intf_supress_4/audio_opt10.raw';
%mic1_file = '../data/invitro_intf_supress_5/audio_opt7.raw';
%mic2_file = '../data/invitro_intf_supress_5/audio_opt11.raw';
%mic1_file = '../data/invitro_intf_supress_6/audio_opt7.raw';
%mic2_file = '../data/invitro_intf_supress_6/audio_opt11.raw';
fid1 = fopen(mic1_file, 'r');
mic1_data = fread(fid1, inf, 'int16');
fclose(fid1);
fid2 = fopen(mic2_file, 'r');
mic2_data = fread(fid2, inf, 'int16');
fclose(fid2);
% User-specified start and end times (in seconds)
start_time = 7;%; 11; %1;%14;%7;% 3; % Change this as needed
%end_time = 16; %13;  % Change this as needed
% Convert times to sample indices
start_sample = round(start_time * Fs);
end_sample = start_sample+length(t)-1;
mic_signals(1, :) = mic1_data(start_sample:end_sample)'-mean(mic1_data(start_sample:end_sample));
mic_signals(2, :) = mic2_data(start_sample:end_sample)'-mean(mic2_data(start_sample:end_sample));
% Normalize signals to prevent clipping during processing
%mic_signals = mic_signals-mean(mic_signals);
mic_signals(1, :) = mic_signals(1, :) / max(abs(mic_signals(1,:)));
mic_signals(2, :) = mic_signals(2, :) / max(abs(mic_signals(2,:)));
c=26;
%ang_n = -90; %cancel lung
ang_n = 89; %cancel heart
dd=0.012;
%dd=0.006;
%dd=0.005;
% test code section

% Apply the Hilbert Transform to Extend Signals to Complex Domain
%analytic_signals = hilbert(mic_signals')'; % NxT

% STFT parameters
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
beamformed_output = beamformed_output / max(abs(beamformed_output(1000:length(beamformed_output)-1000)));

% Plot Signals
figure;
for i = 1:5
    subplot(5, 2, (i-1)*2 + 1);
    switch i
        case 1
%            plot(t, breathing_with_crackles,'c');
            title('Lung Sound (Source)');
        case 2
%            plot(t, heartbeat_signal,'m');
            title('Heartbeat Signal (Source)');
            ylim([-1 1]);
        case 3
            plot(t, mic_signals(1, :),'r');
            ylim([-1,1]);
            title('Microphone 1 (Before Beamforming)');
        case 4
            plot(t, mic_signals(2, :),'b');
            ylim([-1,1]);
            title('Microphone 2 (Before Beamforming)');
        case 5
            %plot(t, beamformed_output);
            plot(t(1000:length(beamformed_output)-1000), beamformed_output(1000:length(beamformed_output)-1000),'g');
            ylim([-1,1]);
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
     %       spectrogram(breathing_with_crackles, window_length, overlap, nfft, Fs, 'yaxis');
            title('Spectrogram - Lung Source');
        case 2
      %      spectrogram(heartbeat_signal, window_length, overlap, nfft, Fs, 'yaxis');
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

    ylim([0 1]); % Adjust frequency range if needed
end
% Play Sounds
play_start= start_time;%2;
play_end =play_start+duration;%7;
play_duration=play_end - play_start;
%disp('heartbeat sound alone...');
%sound(heartbeat_signal(play_start*Fs:play_end*Fs) * volume_scale, Fs);
%audiowrite('../data/outputs/sound/demo3_heart.wav',heartbeat_signal(play_start*Fs:play_end*Fs) * volume_scale, Fs);
%pause(play_duration + 1);

%disp('Lung sound alone...');
%sound(breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
%audiowrite('../data/outputs/sound/demo3_lung.wav',breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
%pause(play_duration + 1);

disp('Sound at Mic 1 (before beamforming)...');
sound(mic_signals(1, :) * volume_scale, Fs);
%audiowrite('../data/outputs/sound/demo3_mic1_input.wav', mic_signals(1, play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

%disp('Sound at Mic 2 (before beamforming)...');
%sound(mic_signals(2, play_start*Fs:play_end*Fs)) * volume_scale, Fs);
%audiowrite('../data/outputs/sound/mic2_input.wav', mic_signals(2, play_start*Fs:play_end*Fs) * volume_scale, Fs);
%pause(play_duration + 1);

disp('Beamformed output (after beamforming)...');
sound(beamformed_output(500:length(beamformed_output)-500)* volume_scale, Fs);
%audiowrite('../data/outputs/sound/demo3_beamformed_output.wav', beamformed_output(play_start*Fs:play_end*Fs) * volume_scale, Fs);
%pause(play_duration + 1);

%disp('Compare with lung sound alone...');
%sound(breathing_with_crackles(play_start*Fs:play_end*Fs) * volume_scale, Fs);
pause(play_duration + 1);

