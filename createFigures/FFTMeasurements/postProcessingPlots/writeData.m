
%%
% =========================================================================
% Hall FFT

AcronymnName = 'Hall';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsHall128x128_allFramesFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsHall_allFramesFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% PETS FFT

AcronymnName = 'PETS';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsPETS_allFramesFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsPETS_allFramesFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% Stuttgart FFT

AcronymnName = 'Stuttgart';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/Stuttgart_144x192_0011_0310FFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/Stuttgart_144x192_0011_0310FFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% VascoDaGamma FFT

AcronymnName = 'VascoDaGamma';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/VascoDaGamaFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/VascoDaGamaFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% baselineHighway FFT

AcronymnName = 'baselineHighway';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsBaselineHighwayFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsBaselineHighwayFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% baselinePETS FFT

AcronymnName = 'baselinePETS';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsBaselinePetsFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsBaselinePetsFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% thermalPark FFT

AcronymnName = 'thermalPark';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsThermalParkFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsThermalParkFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% cameraJitterTraffic FFT

AcronymnName = 'cameraJitterTraffic';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsCameraJitterFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsCameraJitterFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


%%
% =========================================================================
% dynamicBackgroundCanoe FFT

AcronymnName = 'dynamicBackgroundCanoe';

system(['mkdir ', AcronymnName]);

% Measurements ------------------------------------------------------------
load('../results/ResultsDynamicBackgroundCanoeFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsDynamicBackgroundCanoeFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_CSmeasurements_oracle = [AcronymnName, '/FFTCSmeasurementsOracle.dat'];
File_L1L1meas_oracle       = [AcronymnName, '/FFTL1L1measOracle.dat'];
File_L1L1meas_estim        = [AcronymnName, '/FFTL1L1measEstim.dat'];
File_Measurements_L1L1     = [AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn     = [AcronymnName, '/FFTMeasurementsWarn.dat'];

File_EstimL1L1 = [AcronymnName, '/FFTEstimErrorL1L1.dat'];
File_ReconL1L1 = [AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn = [AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_CSmeasurements_oracle, 'w');
fid2   = fopen(File_L1L1meas_oracle      , 'w');
fid3   = fopen(File_L1L1meas_estim       , 'w');
fid4   = fopen(File_Measurements_L1L1    , 'w');
fid5   = fopen(File_Measurements_Warn    , 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));  
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements_L1L1(i)));
    fprintf(fid5,   '%d %d\n', i, round(num_measurements_Warn(i)));
    
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4);


% Relative Errors ---------------------------------------------------------

fid6   = fopen(File_EstimL1L1, 'w');
fid7   = fopen(File_ReconL1L1, 'w');
fid8   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid6, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid7, '%d %3.3f\n', i, max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid8, '%d %3.3f\n', i, max(log10(reconstruction_error_Warn(i)),-9)+9);
end

fclose(fid6); fclose(fid7); fclose(fid8);
% =========================================================================


