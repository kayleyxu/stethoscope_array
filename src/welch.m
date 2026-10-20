% Welch Powe Spectral Density (dB per Hz) Comparison

%clear; close all; clc;

% Load audio files
%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko Recording Sound Export 02_11_2025 11_12 PM.wav"); %white
%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko Recording Sound Export 02_12_2025 12_09 AM white.wav"); 
%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko white chamber_btspeaker.wav");
%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko white chamber_phone.wav");
%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko white chamber_phone_ext.wav");
%[file1, fs1] = audioread("..\data\stethoscope_rec\Ekowhite chamber_phone_ext with cmm.wav");
[file1, fs1] = audioread("..\data\stethoscope_rec\Eko white chamber_btspeaker_with_cmm.wav");

%[file1, fs1] = audioread("..\data\stethoscope_rec\Eko Recording Sound Export 02_11_2025 11_58 PM.wav"); %pink
%[file2, fs2] = audioread("..\data\stethoscope_rec\whitenoise_mems_4000.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\pinknoise_mems_4000.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\untitled.wav");
[file2, fs2] = audioread("..\data\stethoscope_rec\mems_white_noise_in_chamber_btspeaker.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\mems_white_noise_in_chamber_phone.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\mems_white_noise_in_chamber_phone_ext.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\cmm_mems_white_noise_in_chamber_phone_ext.wav");
%[file2, fs2] = audioread("..\data\stethoscope_rec\cmm_mems_white_noise_in_chamber_phone_ext2.wav");
[file3, fs3] = audioread("..\data\stethoscope_rec\cmm_mems_white_noise_in_chamber_btspeaker_2.wav");
start_time = 3;  % Start from 2 seconds
duration_sec = 10; % Take 10 seconds

% Convert time to sample indices
start_sample1 = round(start_time * fs1);
end_sample1 = start_sample1 + round(duration_sec * fs1) - 1;

start_sample2 = round(start_time * fs2);
end_sample2 = start_sample2 + round(duration_sec * fs2) - 1;

start_sample3 = round(start_time * fs3);
end_sample3 = start_sample2 + round(duration_sec * fs3) - 1;

% Ensure we don't exceed signal length
end_sample1 = min(end_sample1, length(file1));
end_sample2 = min(end_sample2, length(file2));
end_sample3 = min(end_sample3, length(file3));

% Extract selected segments
file1 = file1(start_sample1:end_sample1);
file2 = file2(start_sample2:end_sample2);
file3 = file3(start_sample3:end_sample3);


% Optimize window, overlap, and nfft for each sampling rate
window1 = hann(512); % For fs1=4000Hz
noverlap1 = 256; % 50% overlap
nfft1 = 1024; % FFT length

%window2 = hann(4096); % For fs2=48000Hz or 38400Hz
%noverlap2 = 2048; % 50% overlap
%nfft2 = 8192; % FFT length
window2 = hann(4096*4); % For fs2=48000Hz or 38400Hz
noverlap2 = 2048*4; % 50% overlap
nfft2 = 8192*4; % FFT length

window3=window2;
noverlap3=noverlap2;
nfft3=nfft2;

% Compute Welch PSD  (power spectral density)
[pxx1, f1] = pwelch(file1, window1, noverlap1, nfft1, fs1);
[pxx2, f2] = pwelch(file2, window2, noverlap2, nfft2, fs2);
[pxx3, f3] = pwelch(file3, window3, noverlap3, nfft3, fs3);

% Limit frequency range to 10Hz - 4000Hz for fair comparison
freq_min = 20; freq_max = 4000;
idx1 = (f1 >= freq_min & f1 <= freq_max);
idx2 = (f2 >= freq_min & f2 <= freq_max);
idx3 = (f3 >= freq_min & f3 <= freq_max);

f1 = f1(idx1); pxx1 = pxx1(idx1);
f2 = f2(idx2); pxx2 = pxx2(idx2);
f3 = f3(idx3); pxx3 = pxx3(idx3);

% Normalize PSDs within this frequency range
pxx1 = pxx1 / max(pxx1);
pxx2 = pxx2 / max(pxx2);
pxx3 = pxx3 / max(pxx3);

% Plot both PSDs on the same figure
figure;
semilogx(f1, 10*log10(pxx1), 'b', 'LineWidth', 1.5); hold on;
semilogx(f2, 10*log10(pxx2), 'r', 'LineWidth', 1.5); 
semilogx(f3, 10*log10(pxx3), 'g', 'LineWidth', 1.5);

grid on;
xlabel('Frequency (Hz) - Log Scale');
ylabel('Normalized Power/Frequency (dB/Hz)');
title('Power Spectral Density');
%legend({'Littmann Stethoscope', 'MP MEMS Mic', 'CMM MEMS Mic'});

xlim([freq_min,freq_max]);

figure;
plot(f1, 10*log10(pxx1), 'b', 'LineWidth', 1.5); hold on;
plot(f2, 10*log10(pxx2), 'r', 'LineWidth', 1.5);
plot(f3, 10*log10(pxx3), 'g', 'LineWidth', 1.5);
xlim([freq_min,freq_max]);

grid on;
xlabel('Frequency (Hz)');
ylabel('Normalized Power/Frequency (dB/Hz)');
title('Power Spectral Density (Linear Scale)');
legend({'Littmann Stethoscope Recording', 'MP MEMS Mic', 'CMM MEMS Mic'});