% ========================================================================
% DiffBF_sim_and_invitro.m
%
%
% This single script:
%   1. Runs the full-precision Section A simulation (A1 azimuth, A2 depth,
%      A3 common-noise sweep).
%   2. Loads the actual Section B in-vitro recordings.
%   3. Calculates noise-corrected SIR before/after beamforming and DeltaSIR
%      for five conditions with N = 3 repeated acquisitions.
%   4. Creates matched Section A / Section B comparison tables.
%   5. Exports paper-quality PNG, vector PDF, FIG, CSV, LaTeX-table, and MAT.
%
% Condition mapping:
%   cond01 = (-120 deg, 5 cm)
%   cond02 = ( -90 deg, 5 cm)
%   cond03 = ( -60 deg, 5 cm)
%   cond04 = ( -90 deg, 3 cm)
%   cond05 = ( -90 deg, 7 cm)
%
% Required: mic_preprocess, delayseq, diffbfweights, Signal Processing Toolbox
% ========================================================================
clear; clc; close all;
%% ============================= OUTPUT FOLDERS =============================
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
combinedOutputFolder = fullfile( ...
    scriptDir,'results_sectionA_and_actual_invitro');
figureFolder = fullfile(combinedOutputFolder,'figures');
tableFolder = fullfile(combinedOutputFolder,'tables');
if ~exist(figureFolder,'dir'), mkdir(figureFolder); end
if ~exist(tableFolder,'dir'), mkdir(tableFolder); end
%% ============================== DATA FOLDERS ===============================
CRACKLE_FOLDER_OVERRIDE = ...
    '..\data\2024_11_26_23_51_25_single_crackle';
IN_VITRO_WAVE_FOLDER_OVERRIDE = '';
crackleCandidates = { ...
    CRACKLE_FOLDER_OVERRIDE, ...
    fullfile(scriptDir,'..','data','2024_11_26_23_51_25_single_crackle')};
crackleFolderA = '';
for iFolder = 1:numel(crackleCandidates)
    if ~isempty(crackleCandidates{iFolder}) && isfolder(crackleCandidates{iFolder})
        crackleFolderA = crackleCandidates{iFolder}; break;
    end
end
if isempty(crackleFolderA)
    error(['Could not locate the crackle source folder.\n' ...
        'Set CRACKLE_FOLDER_OVERRIDE near the top of this script.']);
end
waveFolderCandidates = { ...
    IN_VITRO_WAVE_FOLDER_OVERRIDE, ...
    fullfile(scriptDir,'..','data','SectionB_invitro_recording','waves'), ...
    fullfile(scriptDir,'..','data','sectionB_invitro_recording','waves')};
inVitroWaveFolder = '';
for iFolder = 1:numel(waveFolderCandidates)
    if ~isempty(waveFolderCandidates{iFolder}) && isfolder(waveFolderCandidates{iFolder})
        inVitroWaveFolder = waveFolderCandidates{iFolder}; break;
    end
end
if isempty(inVitroWaveFolder)
    error(['Could not locate the actual in-vitro WAV folder.\n' ...
        'Set IN_VITRO_WAVE_FOLDER_OVERRIDE near the top of this script.']);
end
fprintf('Crackle folder:\n  %s\n',crackleFolderA);
fprintf('Actual in-vitro WAV folder:\n  %s\n\n',inVitroWaveFolder);
%% ========================================================================
%  SECTION A: FULL-PRECISION SIMULATION
%  ========================================================================
FsA = 44100;
gA = round(0.25*FsA);
rng(0);
[lungA,heartRawA] = A_make_sources(FsA,crackleFolderA);
heartA = heartRawA;
sIdealA = lungA;
[floA,fhiA] = A_lung_band(lungA,FsA,0.95);
analysisBandA = [floA fhiA];
fprintf('lung 95%%-energy band: %.6f - %.6f Hz\n',floA,fhiA);
cA = 40;
ddA = 0.005;
lungPosA = [0.03,0.06,0.05];
heartPosA = [0,-0.07,0.05];
heartWallYA = heartPosA(2);
heartWallSlopeYA = 0;
targetSectionAInputSirDb = 0;
[heartA,heartCalibrationGainA,actualSectionAInputSirDb] = ...
    A_calibrate_heart_to_sir(lungA,heartA,targetSectionAInputSirDb, ...
        lungPosA,heartPosA,cA,ddA,FsA,gA,analysisBandA);
NULLA = A_verify_null_side(lungA,heartA,FsA,cA,ddA, ...
    lungPosA,heartPosA,analysisBandA);
fprintf('empirical null convention: ang_n = %+d nulls the -90 heart side\n',NULLA);
fprintf('Section A heart calibration gain: %.8f\n',heartCalibrationGainA);
fprintf('Section A nominal SIR target/actual: %.4f / %.4f dB\n\n', ...
    targetSectionAInputSirDb,actualSectionAInputSirDb);
fprintf('################## SECTION A: SIMULATION (c=40, dd=0.005) ##################\n');
%% -------------------------- A1 AZIMUTH SWEEP -----------------------------
azlist = [-140 -135 -130 -120 -110 -100 -90 ...
     -80  -70  -60  -50  -45  -40];
A1 = zeros(numel(azlist),7);
A1LungGain = zeros(numel(azlist),1);
A1HeartAtten = zeros(numel(azlist),1);
fprintf('\n-- A1 azimuth sweep (noiseless) --\n');
fprintf('  az |      u       |    SIRin |   SIRout |    dSIR | LungGain | HeartAtten |   SDRout |    dSDR\n');
for i = 1:numel(azlist)
    p = B_heart_position_from_wall(azlist(i),5,heartWallYA,heartWallSlopeYA);
    R = A_metrics(lungA,heartA,sIdealA,lungPosA,p,cA,ddA,NULLA,FsA,gA,[],[],analysisBandA);
    u = A_ueff(p);
    A1(i,:) = [u,R.SIR_in,R.SIR_out,R.SIR_d,R.SDR_in,R.SDR_out,R.SDR_d];
    A1LungGain(i) = R.lung_gain;
    A1HeartAtten(i) = R.heart_atten;
    fprintf(' %4.0f | %+10.6f | %+8.4f | %+8.4f | %+8.4f | %+8.4f | %+10.4f | %+8.4f | %+8.4f\n', ...
        azlist(i),u,R.SIR_in,R.SIR_out,R.SIR_d,R.lung_gain,R.heart_atten,R.SDR_out,R.SDR_d);
end
sectionAAzimuthTable = table(azlist(:),A1(:,1),A1(:,2),A1(:,3),A1(:,4), ...
    A1LungGain,A1HeartAtten,A1(:,5),A1(:,6),A1(:,7), ...
    'VariableNames',{'AzimuthDeg','EffectiveDirectionU', ...
    'SIRin_dB','SIRout_dB','DeltaSIR_dB','LungGain_dB','HeartAttenuation_dB', ...
    'SDRin_dB','SDRout_dB','DeltaSDR_dB'});
%% --------------------------- A2 DEPTH SWEEP ------------------------------
zlist = [2 3 4 5 6 7 8 10];
A2 = zeros(numel(zlist),7);
fprintf('\n-- A2 depth sweep (az=-90, noiseless) --\n');
fprintf(' z(cm)|      u       |    SIRin |   SIRout |    dSIR |    SDRin |   SDRout |    dSDR\n');
for i = 1:numel(zlist)
    p = B_heart_position_from_wall(-90,zlist(i),heartWallYA,heartWallSlopeYA);
    R = A_metrics(lungA,heartA,sIdealA,lungPosA,p,cA,ddA,NULLA,FsA,gA,[],[],analysisBandA);
    u = A_ueff(p);
    A2(i,:) = [u,R.SIR_in,R.SIR_out,R.SIR_d,R.SDR_in,R.SDR_out,R.SDR_d];
    fprintf(' %4.1f | %+10.6f | %+8.4f | %+8.4f | %+8.4f | %+8.4f | %+8.4f | %+8.4f\n', ...
        zlist(i),u,R.SIR_in,R.SIR_out,R.SIR_d,R.SDR_in,R.SDR_out,R.SDR_d);
end
sectionADepthTable = table(zlist(:),A2(:,1),A2(:,2),A2(:,3),A2(:,4), ...
    A2(:,5),A2(:,6),A2(:,7), ...
    'VariableNames',{'DepthCm','EffectiveDirectionU', ...
    'SIRin_dB','SIRout_dB','DeltaSIR_dB','SDRin_dB','SDRout_dB','DeltaSDR_dB'});
%% ------------------- A3 NOISE SWEEP (COMMON NOISE ONLY) ------------------
% Environmental (spatially correlated) noise is the operative case in the
% high-SNR deployed regime. The differential array SUPPRESSES it. Independent
% sensor noise is omitted here and disclosed in text (white-noise-gain cost,
% negligible given the 65 dB MEMS SNR floor).
snrlist = [40 35 30 25 20 15 10];
NrepA = 20;
% Columns: 1-3 = SNRin, SNRout, DeltaSNR ; 4-6 = SINRin, SINRout, DeltaSINR
A3c = zeros(numel(snrlist),6);
A3cS = zeros(numel(snrlist),6);
fprintf('\n-- A3 noise sweep (heart -90), band-limited, partially correlated ambient (alpha=0.8), 10-40 dB, N=%d --\n',NrepA);
fprintf(' inSNR | SNRout(dSNR +/- sd), SINRout(dSINR +/- sd)\n');
for i = 1:numel(snrlist)
    noiseRuns = zeros(NrepA,6);
    for r = 1:NrepA
        Rc = A_metrics(lungA,heartA,sIdealA,lungPosA,heartPosA,cA,ddA,NULLA,FsA,gA, ...
            snrlist(i),[],[floA fhiA],'ambient',1);
        noiseRuns(r,:) = [Rc.SNR_in,Rc.SNR_out,Rc.SNR_d, ...
            Rc.SINR_in,Rc.SINR_out,Rc.SINR_d];
    end
    A3c(i,:) = mean(noiseRuns,1);
    A3cS(i,:) = std(noiseRuns,0,1);
    fprintf(' %5.0f | %8.4f (%+8.4f +/- %.4f), %8.4f (%+8.4f +/- %.4f)\n', ...
        snrlist(i),A3c(i,2),A3c(i,3),A3cS(i,3),A3c(i,5),A3c(i,6),A3cS(i,6));
end
sectionANoiseTable = table(snrlist(:), ...
    A3c(:,1),A3c(:,2),A3c(:,3),A3cS(:,1),A3cS(:,2),A3cS(:,3), ...
    A3c(:,4),A3c(:,5),A3c(:,6),A3cS(:,4),A3cS(:,5),A3cS(:,6), ...
    'VariableNames',{'TargetInputSNR_dB', ...
    'CommonSNRinMean_dB','CommonSNRoutMean_dB','CommonDeltaSNRMean_dB', ...
    'CommonSNRinStd_dB','CommonSNRoutStd_dB','CommonDeltaSNRStd_dB', ...
    'CommonSINRinMean_dB','CommonSINRoutMean_dB','CommonDeltaSINRMean_dB', ...
    'CommonSINRinStd_dB','CommonSINRoutStd_dB','CommonDeltaSINRStd_dB'});
%% -------------------- WRITE SIM CSVs (no supplement figure) --------------
writetable(sectionAAzimuthTable,fullfile(combinedOutputFolder,'sectionA_A1_azimuth_full_precision.csv'));
writetable(sectionADepthTable,fullfile(combinedOutputFolder,'sectionA_A2_depth_full_precision.csv'));
writetable(sectionANoiseTable,fullfile(combinedOutputFolder,'sectionA_A3_noise_common.csv'));

%% ========================================================================
%  SECTION B: ACTUAL IN-VITRO RECORDINGS
%  ========================================================================
nConditionsB = 5;
nRepeatsB = 3;
conditionPairsB = [-120,5; -90,5; -60,5; -90,3; -90,7];
conditionLabelB = ["Off-axis left";"Cardiac-facing center";"Off-axis right";"Shallow";"Deep"];
idxAzLeftB = 1; idxAzCenterB = 2; idxAzRightB = 3;
idxDepthShallowB = 4; idxDepthCenterB = 2; idxDepthDeepB = 5;
cB = 30;                 % measured phantom sound speed (distinct from c=40 thorax model)
ddB = 0.005;
nullAngleB = -89;
bfBandB = [1 3000];
analysisBandB = analysisBandA;
trimSecondsB = 0.25;
fprintf('\n============================================================\n');
fprintf('SECTION B: ACTUAL IN-VITRO RECORDINGS\n');
fprintf('============================================================\n');
fprintf('Analysis band: %.6f - %.6f Hz\n',analysisBandB(1),analysisBandB(2));
fprintf('Beamformer: c = %.3f m/s, d = %.3f mm, null = %.1f deg\n\n',cB,1000*ddB,nullAngleB);
%% ========================== BUILD FILE MANIFEST ===========================
nRowsB = nConditionsB*nRepeatsB;
conditionIndexColumnB = repelem((1:nConditionsB).',nRepeatsB);
repeatColumnB = repmat((1:nRepeatsB).',nConditionsB,1);
lungWavFileB = strings(nRowsB,1);
heartWavFileB = strings(nRowsB,1);
noiseWavFileB = strings(nRowsB,1);
rowB = 0;
for ic = 1:nConditionsB
    for rep = 1:nRepeatsB
        rowB = rowB+1;
        lungWavFileB(rowB) = string(B_find_recording_file(inVitroWaveFolder,ic,rep,'lung'));
        heartWavFileB(rowB) = string(B_find_recording_file(inVitroWaveFolder,ic,rep,'heart'));
        noiseWavFileB(rowB) = string(B_find_recording_file(inVitroWaveFolder,ic,rep,'noise'));
    end
end
sectionBInVitroManifest = table(conditionIndexColumnB, ...
    conditionLabelB(conditionIndexColumnB), ...
    conditionPairsB(conditionIndexColumnB,1), ...
    conditionPairsB(conditionIndexColumnB,2), ...
    repeatColumnB,lungWavFileB,heartWavFileB,noiseWavFileB, ...
    'VariableNames',{'ConditionIndex','ConditionLabel','AzimuthDeg','DepthCm', ...
    'Repeat','LungWavFile','HeartWavFile','NoiseWavFile'});
%% ============================= LOAD + ANALYZE ==============================
sirInRepB = zeros(nConditionsB,nRepeatsB);
sirOutRepB = zeros(nConditionsB,nRepeatsB);
deltaSirRepB = zeros(nConditionsB,nRepeatsB);
floorFlagRepB = false(nConditionsB,nRepeatsB);
sampleRateRepB = zeros(nConditionsB,nRepeatsB);
durationLungSecRepB = zeros(nConditionsB,nRepeatsB);
durationHeartSecRepB = zeros(nConditionsB,nRepeatsB);
durationNoiseSecRepB = zeros(nConditionsB,nRepeatsB);
peakLungRepB = zeros(nConditionsB,nRepeatsB);
peakHeartRepB = zeros(nConditionsB,nRepeatsB);
peakNoiseRepB = zeros(nConditionsB,nRepeatsB);
possibleClippingRepB = false(nConditionsB,nRepeatsB);
referenceFsB = [];
for rowB = 1:nRowsB
    ic = sectionBInVitroManifest.ConditionIndex(rowB);
    rep = sectionBInVitroManifest.Repeat(rowB);
    [lungRead,FsLung] = audioread(char(sectionBInVitroManifest.LungWavFile(rowB)));
    [heartRead,FsHeart] = audioread(char(sectionBInVitroManifest.HeartWavFile(rowB)));
    [noiseRead,FsNoise] = audioread(char(sectionBInVitroManifest.NoiseWavFile(rowB)));
    if FsLung~=FsHeart || FsLung~=FsNoise
        error(['Sample-rate mismatch within condition %d, repeat %d:\n' ...
            'lung=%d Hz, heart=%d Hz, noise=%d Hz.'],ic,rep,FsLung,FsHeart,FsNoise);
    end
    if isempty(referenceFsB)
        referenceFsB = FsLung;
    elseif FsLung~=referenceFsB
        error(['Sample-rate mismatch across recordings. Expected %d Hz but found %d Hz.'], ...
            referenceFsB,FsLung);
    end
    if size(lungRead,2)~=2 || size(heartRead,2)~=2 || size(noiseRead,2)~=2
        error(['Expected stereo recordings for condition %d, repeat %d.\n' ...
            'Found lung=%d, heart=%d, noise=%d channels.'],ic,rep, ...
            size(lungRead,2),size(heartRead,2),size(noiseRead,2));
    end
    FsB = FsLung;
    trimSamplesB = round(trimSecondsB*FsB);
    if analysisBandB(2)>=FsB/2
        error(['The analysis-band upper edge %.3f Hz must be below the WAV Nyquist %.3f Hz.'], ...
            analysisBandB(2),FsB/2);
    end
    if size(lungRead,1)<=2*trimSamplesB || size(heartRead,1)<=2*trimSamplesB || ...
            size(noiseRead,1)<=2*trimSamplesB
        error(['A recording is too short for %.3f s trimming at condition %d, repeat %d.'], ...
            trimSecondsB,ic,rep);
    end
    lungMeasuredMic = lungRead.';
    heartMeasuredMic = heartRead.';
    noiseOnlyMic = noiseRead.';
    [sirResult,floorFlags] = B_calculate_noise_corrected_sir( ...
        lungMeasuredMic,heartMeasuredMic,noiseOnlyMic, ...
        FsB,analysisBandB,trimSamplesB,cB,ddB,nullAngleB,bfBandB);
    sirInRepB(ic,rep) = sirResult(1);
    sirOutRepB(ic,rep) = sirResult(2);
    deltaSirRepB(ic,rep) = sirResult(3);
    floorFlagRepB(ic,rep) = any(floorFlags);
    sampleRateRepB(ic,rep) = FsB;
    durationLungSecRepB(ic,rep) = size(lungRead,1)/FsB;
    durationHeartSecRepB(ic,rep) = size(heartRead,1)/FsB;
    durationNoiseSecRepB(ic,rep) = size(noiseRead,1)/FsB;
    peakLungRepB(ic,rep) = max(abs(lungRead),[],'all');
    peakHeartRepB(ic,rep) = max(abs(heartRead),[],'all');
    peakNoiseRepB(ic,rep) = max(abs(noiseRead),[],'all');
    possibleClippingRepB(ic,rep) = peakLungRepB(ic,rep)>=0.999 || ...
        peakHeartRepB(ic,rep)>=0.999 || peakNoiseRepB(ic,rep)>=0.999;
    fprintf('cond%02d rep%02d: SIRin=%+.4f, SIRout=%+.4f, DeltaSIR=%+.4f dB', ...
        ic,rep,sirResult(1),sirResult(2),sirResult(3));
    if floorFlagRepB(ic,rep), fprintf('  [POWER-SUBTRACTION FLOOR]'); end
    if possibleClippingRepB(ic,rep), fprintf('  [POSSIBLE CLIPPING]'); end
    fprintf('\n');
end
FsB = referenceFsB;
trimSamplesB = round(trimSecondsB*FsB);
%% ============================ IN-VITRO TABLES =============================
sirInMeanB = mean(sirInRepB,2); sirInSDB = std(sirInRepB,0,2);
sirOutMeanB = mean(sirOutRepB,2); sirOutSDB = std(sirOutRepB,0,2);
deltaSirMeanB = mean(deltaSirRepB,2); deltaSirSDB = std(deltaSirRepB,0,2);
anyFloorFlagB = any(floorFlagRepB,2);
anyPossibleClippingB = any(possibleClippingRepB,2);
conditionIndexRawB = repelem((1:nConditionsB).',nRepeatsB);
repeatIndexRawB = repmat((1:nRepeatsB).',nConditionsB,1);
sectionBInVitroRawTable = table(conditionIndexRawB, ...
    conditionLabelB(conditionIndexRawB), ...
    conditionPairsB(conditionIndexRawB,1),conditionPairsB(conditionIndexRawB,2), ...
    repeatIndexRawB,reshape(sampleRateRepB.',[],1), ...
    reshape(durationLungSecRepB.',[],1),reshape(durationHeartSecRepB.',[],1), ...
    reshape(durationNoiseSecRepB.',[],1),reshape(peakLungRepB.',[],1), ...
    reshape(peakHeartRepB.',[],1),reshape(peakNoiseRepB.',[],1), ...
    reshape(sirInRepB.',[],1),reshape(sirOutRepB.',[],1),reshape(deltaSirRepB.',[],1), ...
    reshape(floorFlagRepB.',[],1),reshape(possibleClippingRepB.',[],1), ...
    'VariableNames',{'ConditionIndex','ConditionLabel','AzimuthDeg','DepthCm','Repeat', ...
    'SampleRate_Hz','LungDuration_s','HeartDuration_s','NoiseDuration_s', ...
    'LungPeak','HeartPeak','NoisePeak','SIRin_dB','SIRout_dB','DeltaSIR_dB', ...
    'AnyFloorFlag','PossibleClipping'});
sectionBInVitroSummaryTable = table((1:nConditionsB).', ...
    repmat(nRepeatsB,nConditionsB,1),conditionLabelB, ...
    conditionPairsB(:,1),conditionPairsB(:,2), ...
    sirInMeanB,sirInSDB,sirOutMeanB,sirOutSDB,deltaSirMeanB,deltaSirSDB, ...
    anyFloorFlagB,anyPossibleClippingB, ...
    'VariableNames',{'ConditionIndex','N','ConditionLabel','AzimuthDeg','DepthCm', ...
    'SIRinMean_dB','SIRinSD_dB','SIRoutMean_dB','SIRoutSD_dB', ...
    'DeltaSIRMean_dB','DeltaSIRSD_dB','AnyFloorFlag','PossibleClipping'});
sectionBInVitroAzimuthTable = sectionBInVitroSummaryTable(1:3,:);
sectionBInVitroDepthTable = sectionBInVitroSummaryTable( ...
    [idxDepthShallowB idxDepthCenterB idxDepthDeepB],:);
centerAdvantageDeltaSirB = deltaSirMeanB(idxAzCenterB) - ...
    max(deltaSirMeanB([idxAzLeftB idxAzRightB]));
centerAdvantageSirOutB = sirOutMeanB(idxAzCenterB) - ...
    max(sirOutMeanB([idxAzLeftB idxAzRightB]));
depthTrendCorrectB = deltaSirMeanB(idxDepthShallowB) > deltaSirMeanB(idxDepthCenterB) && ...
    deltaSirMeanB(idxDepthCenterB) > deltaSirMeanB(idxDepthDeepB);
%% ====================== MATCHED SECTION A EXTRACTION ======================
sectionAMatchedSirIn = zeros(nConditionsB,1);
sectionAMatchedSirOut = zeros(nConditionsB,1);
sectionAMatchedDeltaSir = zeros(nConditionsB,1);
for ic = 1:nConditionsB
    azimuthHere = conditionPairsB(ic,1);
    depthHere = conditionPairsB(ic,2);
    if depthHere==5
        idxA = find(sectionAAzimuthTable.AzimuthDeg==azimuthHere,1);
        if isempty(idxA), error('Section A azimuth condition was not found.'); end
        sectionAMatchedSirIn(ic) = sectionAAzimuthTable.SIRin_dB(idxA);
        sectionAMatchedSirOut(ic) = sectionAAzimuthTable.SIRout_dB(idxA);
        sectionAMatchedDeltaSir(ic) = sectionAAzimuthTable.DeltaSIR_dB(idxA);
    else
        idxA = find(sectionADepthTable.DepthCm==depthHere,1);
        if isempty(idxA), error('Section A depth condition was not found.'); end
        sectionAMatchedSirIn(ic) = sectionADepthTable.SIRin_dB(idxA);
        sectionAMatchedSirOut(ic) = sectionADepthTable.SIRout_dB(idxA);
        sectionAMatchedDeltaSir(ic) = sectionADepthTable.DeltaSIR_dB(idxA);
    end
end
bMinusADeltaSir = deltaSirMeanB-sectionAMatchedDeltaSir;
sectionAInVitroComparisonTable = table((1:nConditionsB).',conditionLabelB, ...
    conditionPairsB(:,1),conditionPairsB(:,2), ...
    sectionAMatchedSirIn,sectionAMatchedSirOut,sectionAMatchedDeltaSir, ...
    sirInMeanB,sirInSDB,sirOutMeanB,sirOutSDB,deltaSirMeanB,deltaSirSDB,bMinusADeltaSir, ...
    'VariableNames',{'ConditionIndex','ConditionLabel','AzimuthDeg','DepthCm', ...
    'SectionA_SIRin_dB','SectionA_SIRout_dB','SectionA_DeltaSIR_dB', ...
    'InVitro_SIRinMean_dB','InVitro_SIRinSD_dB','InVitro_SIRoutMean_dB','InVitro_SIRoutSD_dB', ...
    'InVitro_DeltaSIRMean_dB','InVitro_DeltaSIRSD_dB','InVitroMinusA_DeltaSIR_dB'});
centerAdvantageDeltaSirA = sectionAMatchedDeltaSir(idxAzCenterB) - ...
    max(sectionAMatchedDeltaSir([idxAzLeftB idxAzRightB]));
depthTrendCorrectA = sectionAMatchedDeltaSir(idxDepthShallowB) > sectionAMatchedDeltaSir(idxDepthCenterB) && ...
    sectionAMatchedDeltaSir(idxDepthCenterB) > sectionAMatchedDeltaSir(idxDepthDeepB);
overallSummaryTable = table(mean(sectionAMatchedDeltaSir),mean(deltaSirMeanB), ...
    mean(bMinusADeltaSir),centerAdvantageDeltaSirA,centerAdvantageDeltaSirB, ...
    depthTrendCorrectA,depthTrendCorrectB,any(anyFloorFlagB),any(anyPossibleClippingB), ...
    'VariableNames',{'SectionA_MeanDeltaSIR_dB','InVitro_MeanDeltaSIR_dB', ...
    'InVitroMinusA_MeanDeltaSIR_dB','SectionA_CenterAdvantage_dB','InVitro_CenterAdvantage_dB', ...
    'SectionA_DepthTrendCorrect','InVitro_DepthTrendCorrect','AnyInVitroFloorFlag','AnyPossibleClipping'});
fprintf('\n============================================================\n');
fprintf('ACTUAL IN-VITRO MEAN +/- SAMPLE SD, N = %d\n',nRepeatsB);
fprintf('============================================================\n');
disp(sectionBInVitroSummaryTable);
fprintf('\n============================================================\n');
fprintf('MATCHED SECTION A / ACTUAL IN-VITRO COMPARISON\n');
fprintf('============================================================\n');
disp(sectionAInVitroComparisonTable);
fprintf('\n============================================================\n');
fprintf('OVERALL PAPER SUMMARY\n');
fprintf('============================================================\n');
disp(overallSummaryTable);
%% ============================ SAVE ALL TABLES =============================
writetable(sectionAAzimuthTable,fullfile(tableFolder,'SectionA_A1_azimuth.csv'));
writetable(sectionADepthTable,fullfile(tableFolder,'SectionA_A2_depth.csv'));
writetable(sectionANoiseTable,fullfile(tableFolder,'SectionA_A3_noise_common.csv'));
writetable(sectionBInVitroManifest,fullfile(tableFolder,'SectionB_invitro_manifest.csv'));
writetable(sectionBInVitroRawTable,fullfile(tableFolder,'SectionB_invitro_raw_N3.csv'));
writetable(sectionBInVitroSummaryTable,fullfile(tableFolder,'SectionB_invitro_summary_N3.csv'));
writetable(sectionAInVitroComparisonTable,fullfile(tableFolder,'Paper_Table_A_vs_invitro.csv'));
writetable(overallSummaryTable,fullfile(tableFolder,'Paper_Table_overall_summary.csv'));
P_write_latex_tables(sectionAInVitroComparisonTable,sectionBInVitroSummaryTable,tableFolder);
%% ============= PAPER FIGURE 1: SIMULATION (1x3) ==========================
paperFigure1 = figure('Color','w','Units','inches','Position',[1 1 7.2 2.6]);
layout1 = tiledlayout(paperFigure1,1,3,'TileSpacing','compact','Padding','compact');
nexttile(layout1);
plot(azlist,A1(:,4),'-o','LineWidth',1.6,'MarkerSize',5); hold on; xline(-90,'--');
xlabel('heart-source azimuth (deg)'); ylabel('\DeltaSIR (dB)');
title('(a) Azimuth response'); grid on;
nexttile(layout1);
plot(zlist,A2(:,4),'-o','LineWidth',1.6,'MarkerSize',5);
xlabel('heart-source depth (cm)'); ylabel('\DeltaSIR (dB)');
title('(b) Depth response'); grid on;
nexttile(layout1);
errorbar(snrlist,A3c(:,4),A3cS(:,4),'-o','LineWidth',1.4,'MarkerSize',5); hold on;
errorbar(snrlist,A3c(:,5),A3cS(:,5),'-s','LineWidth',1.4,'MarkerSize',5);
xlabel('input lung SNR (dB)'); ylabel('SINR (dB)');
title('(c) Common-noise robustness');
legend('SINR in','SINR out','Location','best'); grid on;
P_style_figure(paperFigure1);
P_export_figure(paperFigure1,figureFolder,'Fig1_SectionA_simulation');

%% ===== PAPER FIGURE 2: SIM-vs-IN-VITRO dSIR OVERLAY (1x2, shared y) ======
yl = [8 12];   % shared y-limits for direct azimuth/depth comparison
paperFigure2 = figure('Color','w','Units','inches','Position',[1 1 7.2 3.0]);
layout2 = tiledlayout(paperFigure2,1,2,'TileSpacing','compact','Padding','compact');
% (a) azimuth: full sim sweep line + in-vitro dots at -120/-90/-60, depth 5 cm
nexttile(layout2);
plot(azlist,A1(:,4),'-','LineWidth',1.6); hold on;
errorbar(sectionBInVitroAzimuthTable.AzimuthDeg, ...
    sectionBInVitroAzimuthTable.DeltaSIRMean_dB, ...
    sectionBInVitroAzimuthTable.DeltaSIRSD_dB, ...
    'o','LineWidth',1.5,'MarkerSize',7);
xline(-90,'--');
xlabel('heart-source azimuth (deg)'); ylabel('\DeltaSIR (dB)');
title('(a) Azimuth, depth 5 cm'); ylim(yl); grid on;
legend('simulation','in-vitro','Location','south');
% (b) depth: full sim sweep line + in-vitro dots at 3/5/7 cm, azimuth -90
nexttile(layout2);
plot(zlist,A2(:,4),'-','LineWidth',1.6); hold on;
errorbar(sectionBInVitroDepthTable.DepthCm, ...
    sectionBInVitroDepthTable.DeltaSIRMean_dB, ...
    sectionBInVitroDepthTable.DeltaSIRSD_dB, ...
    'o','LineWidth',1.5,'MarkerSize',7);
xlabel('heart-source depth (cm)'); ylabel('\DeltaSIR (dB)');
title('(b) Depth, azimuth -90 deg'); ylim(yl); grid on;
legend('simulation','in-vitro','Location','south');
P_style_figure(paperFigure2);
P_export_figure(paperFigure2,figureFolder,'Fig2_sim_vs_invitro_overlay');

%% ================================ SAVE MAT ================================
save(fullfile(combinedOutputFolder,'sectionA_and_actual_invitro_all_paper_results.mat'), ...
    'A1','A2','A3c','A3cS', ...
    'sectionAAzimuthTable','sectionADepthTable','sectionANoiseTable', ...
    'sectionBInVitroManifest','sectionBInVitroRawTable','sectionBInVitroSummaryTable', ...
    'sectionBInVitroAzimuthTable','sectionBInVitroDepthTable', ...
    'sectionAInVitroComparisonTable','overallSummaryTable', ...
    'sirInRepB','sirOutRepB','deltaSirRepB','floorFlagRepB','possibleClippingRepB', ...
    'sectionAMatchedSirIn','sectionAMatchedSirOut','sectionAMatchedDeltaSir','bMinusADeltaSir', ...
    'centerAdvantageDeltaSirA','centerAdvantageDeltaSirB','centerAdvantageSirOutB', ...
    'depthTrendCorrectA','depthTrendCorrectB','analysisBandA','analysisBandB', ...
    'cA','cB','ddA','ddB','NULLA','nullAngleB','conditionPairsB','conditionLabelB', ...
    'inVitroWaveFolder','crackleFolderA');
fprintf('\n============================================================\n');
fprintf('PAPER RESULTS COMPLETE\n');
fprintf('============================================================\n');
fprintf('All outputs saved in:\n  %s\n',combinedOutputFolder);
%% =========================== LOCAL FUNCTIONS ==============================
function NULL = A_verify_null_side(lung,heart,Fs,c,dd,lungPos,heartPos,band)
    g = round(0.25*Fs); best = -inf; NULL = -90;
    for candidate = [-90 -89]
        R = A_metrics(lung,heart,lung,lungPos,heartPos,c,dd,candidate,Fs,g,[],[],band);
        if isfinite(R.SIR_out) && R.SIR_out>best, best = R.SIR_out; NULL = candidate; end
    end
end
function [heartScaled,scale,actualSirDb] = A_calibrate_heart_to_sir( ...
        lung,heart,targetSirDb,lungPos,heartPos,c,dd,Fs,g,band)
    micPositions = A_micpos(c,dd,-89,2);
    lungMic = A_buildpos(lung,lungPos,micPositions,c,Fs,[]);
    heartMic = A_buildpos(heart,heartPos,micPositions,c,Fs,[]);
    lungPower = A_bandpower_seg(lungMic(1,:),Fs,band,g);
    heartPower = A_bandpower_seg(heartMic(1,:),Fs,band,g);
    scale = sqrt(lungPower/(max(heartPower,eps)*10^(targetSirDb/10)));
    heartScaled = scale*heart;
    heartMicScaled = A_buildpos(heartScaled,heartPos,micPositions,c,Fs,[]);
    heartPowerScaled = A_bandpower_seg(heartMicScaled(1,:),Fs,band,g);
    actualSirDb = 10*log10(lungPower/max(heartPowerScaled,eps));
end
function R = A_metrics(lung,heart,idealLung,lungPos,heartPos, ...
        c,dd,nullAngle,Fs,g,inputSnrDb,effect,band,noiseType,nNoiseDraws)
    if nargin<13 || isempty(band), band = [1 Fs/2]; end
    if nargin<14 || isempty(noiseType), noiseType = 'indep'; end
    if nargin<15 || isempty(nNoiseDraws), nNoiseDraws = 40; end
    micPositions = A_micpos(c,dd,nullAngle,2);
    lungMic = A_buildpos(lung,lungPos,micPositions,c,Fs,effect);
    heartMic = A_buildpos(heart,heartPos,micPositions,c,Fs,effect);
    lungBf = A_bf(lungMic,Fs,c,dd,nullAngle);
    heartBf = A_bf(heartMic,Fs,c,dd,nullAngle);
    powerBand = @(x) A_bandpower_seg(x,Fs,band,g);
    R.SIR_in = 10*log10(powerBand(lungMic(1,:))/powerBand(heartMic(1,:)));
    R.SIR_out = 10*log10(powerBand(lungBf)/powerBand(heartBf));
    R.SIR_d = R.SIR_out-R.SIR_in;
    R.lung_gain = 10*log10(powerBand(lungBf)/powerBand(lungMic(1,:)));
    R.heart_atten = 10*log10(powerBand(heartMic(1,:))/powerBand(heartBf));
    if ~isempty(inputSnrDb)
        lungBandPower = powerBand(lungMic(1,:));
        noiseInputAccumulator = 0; noiseOutputAccumulator = 0;
        alphaCorr = 0.8;                  % correlated (shared ambient) fraction
        gainMismatch = 10^(0.5/20);       % <=0.5 dB inter-channel gain mismatch
        delaySamples = 8e-6*Fs;           % <=8 us inter-channel delay
        for q = 1:nNoiseDraws
            % Partially correlated ambient field: shared component + per-channel
            % independent component, with small gain/delay mismatch. noiseType is
            % retained for interface compatibility but no longer branches.
            sharedNoise = randn(1,numel(lung));
            indepNoise1 = randn(1,numel(lung));
            indepNoise2 = randn(1,numel(lung));
            sharedDelayed = delayseq(sharedNoise.',delaySamples).';
            unitNoise = [ sqrt(alphaCorr)*sharedNoise + sqrt(1-alphaCorr)*indepNoise1; ...
                          sqrt(alphaCorr)*gainMismatch*sharedDelayed + sqrt(1-alphaCorr)*indepNoise2 ];
            scale = sqrt(lungBandPower/max(powerBand(unitNoise(1,:)),eps)/(10^(inputSnrDb/10)));
            unitNoise = scale*unitNoise;
            noiseBf = A_bf(unitNoise,Fs,c,dd,nullAngle);
            noiseInputAccumulator = noiseInputAccumulator + powerBand(unitNoise(1,:));
            noiseOutputAccumulator = noiseOutputAccumulator + powerBand(noiseBf);
        end
        noiseInputPower = noiseInputAccumulator/nNoiseDraws;
        noiseOutputPower = noiseOutputAccumulator/nNoiseDraws;
        R.SNR_in = 10*log10(powerBand(lungMic(1,:))/max(noiseInputPower,eps));
        R.SNR_out = 10*log10(powerBand(lungBf)/max(noiseOutputPower,eps));
        R.SNR_d = R.SNR_out-R.SNR_in;
        R.SINR_in = 10*log10(powerBand(lungMic(1,:))/max(powerBand(heartMic(1,:))+noiseInputPower,eps));
        R.SINR_out = 10*log10(powerBand(lungBf)/max(powerBand(heartBf)+noiseOutputPower,eps));
        R.SINR_d = R.SINR_out-R.SINR_in;
    else
        R.SNR_in = NaN; R.SNR_out = NaN; R.SNR_d = NaN;
        R.SINR_in = 10*log10(powerBand(lungMic(1,:))/max(powerBand(heartMic(1,:)),eps));
        R.SINR_out = 10*log10(powerBand(lungBf)/max(powerBand(heartBf),eps));
        R.SINR_d = R.SINR_out-R.SINR_in;
    end
    cleanOutputMixture = lungBf+heartBf;
    cleanInputMixture = lungMic(1,:)+heartMic(1,:);
    R.SDR_in = A_sdr_fit(cleanInputMixture,idealLung,g,Fs);
    R.SDR_out = A_sdr_fit(cleanOutputMixture,idealLung,g,Fs);
    R.SDR_d = R.SDR_out-R.SDR_in;
end
function [lowerFrequency,upperFrequency] = A_lung_band(x,Fs,fraction)
    x = x(:)-mean(x); nfft = 2^nextpow2(numel(x));
    spectrumPower = abs(fft(x,nfft)).^2; halfLength = nfft/2;
    spectrumPower = spectrumPower(1:halfLength);
    frequencies = (0:halfLength-1)*Fs/nfft;
    cumulativeFraction = cumsum(spectrumPower)/sum(spectrumPower);
    lowerFraction = (1-fraction)/2; upperFraction = 1-lowerFraction;
    lowerFrequency = frequencies(find(cumulativeFraction>=lowerFraction,1));
    upperFrequency = frequencies(find(cumulativeFraction>=upperFraction,1));
    if isempty(lowerFrequency) || lowerFrequency<20, lowerFrequency = 50; end
    if isempty(upperFrequency), upperFrequency = Fs/2; end
end
function p = A_bandpower_seg(x,Fs,band,g)
    x = x(:); firstSample = g+1; lastSample = numel(x)-g;
    if lastSample<=firstSample, firstSample = 1; lastSample = numel(x); end
    x = x(firstSample:lastSample);
    bp = designfilt('bandpassiir','FilterOrder',8, ...
        'HalfPowerFrequency1',max(1,band(1)),'HalfPowerFrequency2',min(band(2),Fs/2-1), ...
        'SampleRate',Fs);
    y = filtfilt(bp,x); p = sum(y.^2)/numel(y);
end
function micPositions = A_micpos(c,dd,nullAngle,nMics)
    designWavelength = c/100; normalizedSpacing = dd/designWavelength;
    [~,positions] = diffbfweights(nMics,normalizedSpacing,nullAngle,'ArrayGeometry','ULA');
    micPositions = positions'*designWavelength;
end
function u = A_ueff(position)
    radiusXY = hypot(position(1),position(2));
    elevationDeg = atand(position(3)/max(radiusXY,eps));
    azimuthDeg = atan2d(position(2),position(1));
    u = sind(azimuthDeg)*cosd(elevationDeg);
end
function microphoneSignals = A_buildpos(source,sourcePosition,micPositions,c,Fs,effect)
    nMics = size(micPositions,1); microphoneSignals = zeros(nMics,numel(source));
    source = source(:);
    for micIndex = 1:nMics
        propagationDelaySec = norm(micPositions(micIndex,:)-sourcePosition)/c;
        microphoneSignals(micIndex,:) = delayseq(source,propagationDelaySec*Fs).';
    end
    if isequal(effect,'reverb')
        positionKey = mod(round(sum(abs(sourcePosition(:)).*[1;7;13])*1e4),2^31-1);
        rs = RandStream('mt19937ar','Seed',positionKey);
        baseDelaySec = 0.002 + (rand(rs)-0.5)*0.002;
        reflectionLevelDb = -12 + (rand(rs)-0.5)*6;
        reflectionGain = 10^(reflectionLevelDb/20);
        mic1ReflectionDelaySamples = (baseDelaySec + (rand(rs)-0.5)*0.0002)*Fs;
        mic2ReflectionDelaySamples = (baseDelaySec + (rand(rs)-0.5)*0.0002)*Fs;
        microphoneSignals = microphoneSignals + reflectionGain*[ ...
            delayseq(microphoneSignals(1,:).',mic1ReflectionDelaySamples).'; ...
            delayseq(microphoneSignals(2,:).',mic2ReflectionDelaySamples).'];
    end
end
function y = A_bf(microphoneSignals,Fs,c,dd,nullAngle)
    windowLength = 4096; overlapLength = 2048; nfft = 4096;
    minimumFrequency = 1; maximumFrequency = 3000;
    [S1,frequency,~] = stft(microphoneSignals(1,:),Fs, ...
        'Window',hamming(windowLength,'periodic'),'OverlapLength',overlapLength, ...
        'FFTLength',nfft,'Centered',false);
    [S2,~,~] = stft(microphoneSignals(2,:),Fs, ...
        'Window',hamming(windowLength,'periodic'),'OverlapLength',overlapLength, ...
        'FFTLength',nfft,'Centered',false);
    outputSpectrum = zeros(size(S1)); maximumFrequency = min(maximumFrequency,Fs/2);
    for k = 1:length(frequency)
        if frequency(k)>=minimumFrequency && frequency(k)<=maximumFrequency
            wavelength = c/frequency(k); normalizedSpacing = dd/wavelength;
            W = diffbfweights(2,normalizedSpacing,nullAngle,'ArrayGeometry','ULA');
            if any(~isfinite(W))
                outputSpectrum(k,:) = 0.1*(S1(k,:)+S2(k,:));
            else
                outputSpectrum(k,:) = S1(k,:)*W(1)' + S2(k,:)*W(2)';
            end
        else
            outputSpectrum(k,:) = 0.1*(S1(k,:)+S2(k,:));
        end
    end
    outputSpectrum(~isfinite(outputSpectrum)) = 0;
    y = real(istft(outputSpectrum,Fs,'Window',hamming(windowLength,'periodic'), ...
        'OverlapLength',overlapLength,'FFTLength',nfft,'Centered',false));
    y = y(:).';
end
function sdr = A_sdr_fit(x,reference,g,Fs)
    x = x(:); reference = reference(:);
    nCommon = min(numel(x),numel(reference)); x = x(1:nCommon); reference = reference(1:nCommon);
    [correlation,lags] = xcorr(x,reference,round(0.02*Fs));
    [~,maximumIndex] = max(abs(correlation)); lag = lags(maximumIndex);
    alignedReference = zeros(nCommon,1);
    sourceStart = max(1,1-lag); destinationStart = max(1,1+lag);
    nCopy = min(numel(reference)-sourceStart+1,nCommon-destinationStart+1);
    if nCopy>0
        alignedReference(destinationStart:destinationStart+nCopy-1) = ...
            reference(sourceStart:sourceStart+nCopy-1);
    end
    firstSample = g+1; lastSample = nCommon-g;
    if lastSample<=firstSample, firstSample = 1; lastSample = nCommon; end
    xRegion = x(firstSample:lastSample);
    referenceRegion = alignedReference(firstSample:lastSample);
    fittedGain = (referenceRegion'*xRegion)/(referenceRegion'*referenceRegion+eps);
    fittedSignal = fittedGain*referenceRegion; residual = xRegion-fittedSignal;
    sdr = 10*log10(sum(fittedSignal.^2)/(sum(residual.^2)+eps));
end
function [lung,heart] = A_make_sources(Fs,crackleFolder)
    rng(1); durationSec = 10; time = 0:1/Fs:durationSec;
    [crackleArray,~,~] = mic_preprocess(crackleFolder); crackle = crackleArray{1};
    breathingRate = 2/durationSec; cycleDuration = 1/breathingRate;
    inhaleRatio = 0.4; inhaleDuration = inhaleRatio*cycleDuration;
    exhaleDuration = (1-inhaleRatio)*cycleDuration;
    [bLow,aLow] = butter(4,500/(Fs/2),'low'); lung = zeros(size(time));
    for cycle = 1:2
        cycleStart = (cycle-1)*cycleDuration; inhaleEnd = cycleStart+inhaleDuration;
        inhaleMask = (time>=cycleStart) & (time<inhaleEnd); nInhale = sum(inhaleMask);
        inhaleTime = linspace(0,inhaleDuration,nInhale);
        inhaleWave = 0.4*(1-cos(2*pi*inhaleTime/inhaleDuration)).*randn(1,nInhale);
        lung(inhaleMask) = filter(bLow,aLow,inhaleWave);
        exhaleStart = inhaleEnd; exhaleEnd = cycleStart+cycleDuration;
        exhaleMask = (time>=exhaleStart) & (time<exhaleEnd); nExhale = sum(exhaleMask);
        exhaleTime = linspace(0,exhaleDuration,nExhale);
        exhaleWave = 0.8*exp(-exhaleTime/(exhaleDuration/3)).*randn(1,nExhale);
        fadeLength = round(0.1*nExhale); fadeWindow = hann(2*fadeLength).';
        fadeWindow = fadeWindow(1:fadeLength);
        exhaleWave(1:fadeLength) = exhaleWave(1:fadeLength).*fadeWindow;
        lung(exhaleMask) = filter(bLow,aLow,exhaleWave);
        crackleTime = exhaleStart+exhaleDuration/2; crackleStart = round(crackleTime*Fs);
        crackleEnd = crackleStart+numel(crackle)-1;
        if crackleEnd<=numel(time)
            addedCrackle = 15*crackle(1:(crackleEnd-crackleStart+1));
            lung(crackleStart:crackleEnd) = lung(crackleStart:crackleEnd) + addedCrackle(:).';
        end
    end
    lung = lung/max(abs(lung));
    heartRate = 80/60; s1Duration = 0.10; s2Duration = 0.07; systoleDuration = 0.30;
    diastoleDuration = (1/heartRate) - (s1Duration+s2Duration+systoleDuration);
    heartPattern = [ones(1,round(s1Duration*Fs)),zeros(1,round(systoleDuration*Fs)), ...
        ones(1,round(s2Duration*Fs)),zeros(1,round(diastoleDuration*Fs))];
    heart = repmat(heartPattern,1,ceil(numel(time)/numel(heartPattern)));
    heart = circshift(heart(1:numel(time)),round(0.1*Fs));
    heart = bandpass(heart,[50 200],Fs); heart = 0.5*heart/max(abs(heart));
end
function pos = B_heart_position_from_wall(azimuthDeg,depthCm,heartWallY,wallSlopeYperX)
    denominator = tand(azimuthDeg)-wallSlopeYperX;
    if abs(cosd(azimuthDeg))<1e-12, x = 0; else, x = heartWallY/denominator; end
    y = heartWallY + wallSlopeYperX*x;
    pos = [x,y,depthCm/100];
end
function filePath = B_find_recording_file(waveFolder,conditionIndex,repeatIndex,signalType)
    baseName = sprintf('cond%02d_rep%02d_%s',conditionIndex,repeatIndex,signalType);
    directCandidates = {fullfile(waveFolder,[baseName '.wav']),fullfile(waveFolder,[baseName '.WAV'])};
    for i = 1:numel(directCandidates)
        if isfile(directCandidates{i}), filePath = directCandidates{i}; return; end
    end
    listing = dir(waveFolder); targetName = lower([baseName '.wav']);
    for i = 1:numel(listing)
        if ~listing(i).isdir && strcmpi(listing(i).name,targetName)
            filePath = fullfile(waveFolder,listing(i).name); return;
        end
    end
    error('Missing recording: %s.wav in folder:\n%s',baseName,waveFolder);
end
function [result,floorFlags] = B_calculate_noise_corrected_sir( ...
        lungMeasuredMic,heartMeasuredMic,noiseOnlyMic, ...
        Fs,analysisBand,trimSamples,cDesign,ddDesign,nullAngle,bfBand)
    lungMeasuredInputPower = B_band_power_trimmed(lungMeasuredMic(1,:),Fs,analysisBand,trimSamples);
    heartMeasuredInputPower = B_band_power_trimmed(heartMeasuredMic(1,:),Fs,analysisBand,trimSamples);
    noiseInputPower = B_band_power_trimmed(noiseOnlyMic(1,:),Fs,analysisBand,trimSamples);
    lungInputDifference = lungMeasuredInputPower-noiseInputPower;
    heartInputDifference = heartMeasuredInputPower-noiseInputPower;
    floorFlags(1) = lungInputDifference<=0; floorFlags(2) = heartInputDifference<=0;
    correctedLungInputPower = max(lungInputDifference,eps);
    correctedHeartInputPower = max(heartInputDifference,eps);
    sirIn = 10*log10(correctedLungInputPower/correctedHeartInputPower);
    lungBf = B_apply_beamformer(lungMeasuredMic,Fs,cDesign,ddDesign,nullAngle,bfBand);
    heartBf = B_apply_beamformer(heartMeasuredMic,Fs,cDesign,ddDesign,nullAngle,bfBand);
    noiseBf = B_apply_beamformer(noiseOnlyMic,Fs,cDesign,ddDesign,nullAngle,bfBand);
    lungMeasuredOutputPower = B_band_power_trimmed(lungBf,Fs,analysisBand,trimSamples);
    heartMeasuredOutputPower = B_band_power_trimmed(heartBf,Fs,analysisBand,trimSamples);
    noiseOutputPower = B_band_power_trimmed(noiseBf,Fs,analysisBand,trimSamples);
    lungOutputDifference = lungMeasuredOutputPower-noiseOutputPower;
    heartOutputDifference = heartMeasuredOutputPower-noiseOutputPower;
    floorFlags(3) = lungOutputDifference<=0; floorFlags(4) = heartOutputDifference<=0;
    correctedLungOutputPower = max(lungOutputDifference,eps);
    correctedHeartOutputPower = max(heartOutputDifference,eps);
    sirOut = 10*log10(correctedLungOutputPower/correctedHeartOutputPower);
    result = [sirIn,sirOut,sirOut-sirIn];
end
function y = B_apply_beamformer(mic,Fs,cDesign,ddDesign,nullAngle,processingBand)
    win = 4096; overlap = 2048; nfft = 4096; window = hamming(win,'periodic');
    [S1,f,~] = stft(mic(1,:),Fs,'Window',window,'OverlapLength',overlap,'FFTLength',nfft,'Centered',false);
    [S2,~,~] = stft(mic(2,:),Fs,'Window',window,'OverlapLength',overlap,'FFTLength',nfft,'Centered',false);
    Sout = zeros(size(S1));
    for k = 1:numel(f)
        if f(k)>=processingBand(1) && f(k)<=processingBand(2)
            wavelength = cDesign/f(k); normalizedSpacing = ddDesign/wavelength;
            W = diffbfweights(2,normalizedSpacing,nullAngle,'ArrayGeometry','ULA');
            if any(~isfinite(W)), error('Nonfinite beamformer weight at %.2f Hz.',f(k)); end
            Sout(k,:) = S1(k,:)*conj(W(1)) + S2(k,:)*conj(W(2));
        else
            Sout(k,:) = 0.1*(S1(k,:)+S2(k,:));
        end
    end
    y = real(istft(Sout,Fs,'Window',window,'OverlapLength',overlap,'FFTLength',nfft,'Centered',false));
    y = y(:).';
end
function p = B_band_power_trimmed(x,Fs,band,trimSamples)
    x = x(:); firstSample = trimSamples+1; lastSample = numel(x)-trimSamples;
    if lastSample<=firstSample, firstSample = 1; lastSample = numel(x); end
    x = x(firstSample:lastSample); x = x-mean(x);
    bp = designfilt('bandpassiir','FilterOrder',8, ...
        'HalfPowerFrequency1',max(1,band(1)),'HalfPowerFrequency2',min(band(2),Fs/2-1),'SampleRate',Fs);
    y = filtfilt(bp,x); p = mean(y.^2);
end
function P_style_figure(fig)
    axesList = findall(fig,'Type','axes');
    for i = 1:numel(axesList)
        ax = axesList(i);
        set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.8,'Box','off','TickDir','out');
        ax.XGrid = 'on'; ax.YGrid = 'on'; ax.GridAlpha = 0.15; ax.MinorGridAlpha = 0.08;
    end
end
function P_export_figure(fig,figureFolder,baseName)
    figPath = fullfile(figureFolder,[baseName '.fig']);
    pngPath = fullfile(figureFolder,[baseName '.png']);
    pdfPath = fullfile(figureFolder,[baseName '.pdf']);
    savefig(fig,figPath);
    try
        exportgraphics(fig,pngPath,'Resolution',300);
        exportgraphics(fig,pdfPath,'ContentType','vector');
    catch
        print(fig,pngPath,'-dpng','-r300'); print(fig,pdfPath,'-dpdf','-painters');
    end
end
function P_write_latex_tables(comparisonTable,inVitroSummaryTable,tableFolder)
    comparisonPath = fullfile(tableFolder,'Paper_Table_A_vs_invitro.tex');
    fid = fopen(comparisonPath,'w');
    if fid<0, error('Could not create %s.',comparisonPath); end
    cleanupComparison = onCleanup(@() fclose(fid));
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Matched simulation and in-vitro $\\Delta$SIR. In-vitro values are mean $\\pm$ sample SD, $n=3$.}\n');
    fprintf(fid,'\\label{tab:sim_invitro_comparison}\n\\begin{tabular}{lccc}\n\\hline\n');
    fprintf(fid,'Configuration & Simulation & In-vitro & Difference \\\\\n & (dB) & (dB) & (dB) \\\\\n\\hline\n');
    for i = 1:height(comparisonTable)
        fprintf(fid,'$%+.0f^{\\circ}$, %.0f\\,cm & %.2f & $%.2f \\pm %.2f$ & %+.2f \\\\\n', ...
            comparisonTable.AzimuthDeg(i),comparisonTable.DepthCm(i), ...
            comparisonTable.SectionA_DeltaSIR_dB(i),comparisonTable.InVitro_DeltaSIRMean_dB(i), ...
            comparisonTable.InVitro_DeltaSIRSD_dB(i),comparisonTable.InVitroMinusA_DeltaSIR_dB(i));
    end
    fprintf(fid,'\\hline\n\\end{tabular}\n\\end{table}\n');
    clear cleanupComparison;
    inVitroPath = fullfile(tableFolder,'Paper_Table_invitro_before_after.tex');
    fid = fopen(inVitroPath,'w');
    if fid<0, error('Could not create %s.',inVitroPath); end
    cleanupInVitro = onCleanup(@() fclose(fid));
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{In-vitro SIR before and after beamforming. Values are mean $\\pm$ sample SD, $n=3$.}\n');
    fprintf(fid,'\\label{tab:invitro_sir}\n\\begin{tabular}{lccc}\n\\hline\n');
    fprintf(fid,'Configuration & SIR in & SIR out & $\\Delta$SIR \\\\\n & (dB) & (dB) & (dB) \\\\\n\\hline\n');
    for i = 1:height(inVitroSummaryTable)
        fprintf(fid,'$%+.0f^{\\circ}$, %.0f\\,cm & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ \\\\\n', ...
            inVitroSummaryTable.AzimuthDeg(i),inVitroSummaryTable.DepthCm(i), ...
            inVitroSummaryTable.SIRinMean_dB(i),inVitroSummaryTable.SIRinSD_dB(i), ...
            inVitroSummaryTable.SIRoutMean_dB(i),inVitroSummaryTable.SIRoutSD_dB(i), ...
            inVitroSummaryTable.DeltaSIRMean_dB(i),inVitroSummaryTable.DeltaSIRSD_dB(i));
    end
    fprintf(fid,'\\hline\n\\end{tabular}\n\\end{table}\n');
    clear cleanupInVitro;
end


function [cropped_audio_data_20ms, common_window_time_ms, latencies_ms] = mic_preprocess(folder_path, start_time, end_time)
% MIC_PREPROCESS Processes microphone audio files to find and crop the most active 20ms window.
% 
% Inputs:
%   folder_path  - Path to the folder containing the 16 mic audio files in .raw format.
%   start_time   - (Optional) Start time of the first-order cropping in seconds.
%   end_time     - (Optional) End time of the first-order cropping in seconds.
% 
% Outputs:
%   cropped_audio_data_20ms - A cell array (1x16) containing the cropped 20ms audio data for all microphones.
%   common_window_time_ms   - The start time of the 20ms window in milliseconds, relative to the original file.
%   latencies_ms            - Latencies for each mic relative to the earliest mic in ms.
%
% Figures:
%   - Waveforms and delays for all mics (50ms window).
%   - STFT spectrogram for each microphone (50ms window).

%% Constants
num_mics = 5;
Fs = 64000;            % Sampling rate in Hz
bits_per_sample = 16;  % 16-bit audio
sample_format = 'int16'; % Format of raw data (signed 16-bit integers)
% these calibration values are obtained by running source at center scenario
cal_samples=[0     0     0     0     0     0     1     1     2     1     1     1     1     1     0     0];
cal_ms=[0         0         0         0         0         0    0.0208    0.0208    0.0417    0.0208    0.0208    0.0208    0.0208    0.0208         0         0];

%% Step 1: Load and first-order crop the 16 audio files
%disp('Loading and optionally cropping audio files...');
audio_files = cell(num_mics, 1);

% Crop indices (if start_time and end_time are provided)
if nargin > 1 && ~isempty(start_time) && ~isempty(end_time)
    crop_start_idx = round(start_time * Fs);
    crop_end_idx = round(end_time * Fs);
else
    crop_start_idx = 1;
    crop_end_idx = Inf;
end

% Load and optionally crop each file
for mic = 1:num_mics
    
    %temp 3 good mics only
    %{
    if (mic==1) 
        file_name = sprintf('audio_opt0.raw');
    elseif (mic==2)
        file_name = sprintf('audio_opt3.raw');
    else 
        file_name = sprintf('audio_opt4.raw');
    end 
    %}

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
%temp
    
    audio_files{2} = audio_files{1};
    audio_files{3} = audio_files{4};
    


%{
    if nargin > 1
    disp(['First-order cropping applied: Start Time = ', num2str(start_time), ' s, End Time = ', num2str(end_time), ' s']);
else
    disp('No first-order cropping applied (using full signals).');
end
%}
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
%figure;
%sgtitle('STFT Spectrogram for Each Microphone (50ms Window)');
%for mic = 1:num_mics
    % Perform STFT on the 50ms window signal
%    signal = cropped_audio_files_50ms{mic};
%    [S, F, T] = stft(signal, Fs, 'Window', hamming(window_size), ...
%                     'OverlapLength', overlap_size, 'FFTLength', fft_size);
    
    % Magnitude of the STFT
%    magnitude = abs(S);
    
    % Plot the spectrogram
%    subplot(4, 4, mic);
 %   surf(T * 1000, F, 20 * log10(magnitude), 'EdgeColor', 'none'); % Plot in ms
 %   axis tight;
  %  view(0, 90);
  %  colormap jet;
  %  colorbar;
  %  title(['Mic ', num2str(mic-1)]);
  %  xlabel('Time (ms)');
  %  ylabel('Frequency (Hz)');
  %  ylim([0, Fs/2]); % Show only positive frequencies
end

