MotherDir = 'BoxPlotData/';
system(['mkdir ', MotherDir]);

%%
% =========================================================================
% Hall FFT

AcronymnName = 'Hall';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
load('../results/ResultsHall128x128_allFramesFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsHall_allFramesFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% PETS FFT

AcronymnName = 'PETS';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
load('../results/ResultsPETS_allFramesFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsPETS_allFramesFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% Stuttgart FFT

AcronymnName = 'Stuttgart';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
load('../results/Stuttgart_144x192_0011_0310FFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/Stuttgart_144x192_0011_0310FFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% VascoDaGama FFT

AcronymnName = 'VascoDaGama';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
load('../results/VascoDaGamaFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/VascoDaGamaFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% baselineHighway FFT

AcronymnName = 'baselineHighway';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
load('../results/ResultsBaselineHighwayFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsBaselineHighwayFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

% Filenames to save data
File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% baselinePETS FFT

AcronymnName = 'baselinePETS';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
% Filenames to save data
load('../results/ResultsBaselinePetsFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsBaselinePetsFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% thermalPark FFT

AcronymnName = 'thermalPark';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
% Filenames to save data
load('../results/ResultsThermalParkFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsThermalParkFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% cameraJitterTraffic FFT

AcronymnName = 'cameraJitterTraffic';

system(['mkdir ', MotherDir, AcronymnName]);


% Measurements ------------------------------------------------------------
% Filenames to save data
load('../results/ResultsCameraJitterFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsCameraJitterFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%16.16f\n', max(log10(reconstruction_error_L1L1(i)),-18)+18);
    fprintf(fid4, '%16.16f\n', max(log10(reconstruction_error_Warn(i)),-18)+18);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


%%
% =========================================================================
% dynamicBackgroundCanoe FFT

AcronymnName = 'dynamicBackgroundCanoe';

system(['mkdir ', MotherDir, AcronymnName]);

%%
% Measurements ------------------------------------------------------------
% Filenames to save data
load('../results/ResultsDynamicBackgroundCanoeFFT.mat');
num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;
load('../results/ResultsDynamicBackgroundCanoeFFT_Warnell.mat');
num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

File_Measurements_L1L1 = [MotherDir, AcronymnName, '/FFTMeasurementsL1L1.dat'];
File_Measurements_Warn = [MotherDir, AcronymnName, '/FFTMeasurementsWarn.dat'];
File_ReconL1L1         = [MotherDir, AcronymnName, '/FFTReconErrorL1L1.dat'];
File_ReconWarn         = [MotherDir, AcronymnName, '/FFTReconErrorWarn.dat'];

nFrames = Execution_parameters.nFrames;

fid1   = fopen(File_Measurements_L1L1, 'w');
fid2   = fopen(File_Measurements_Warn, 'w');
fid3   = fopen(File_ReconL1L1, 'w');
fid4   = fopen(File_ReconWarn, 'w');

for i = 1 : nFrames - 1 
    fprintf(fid1, '%d\n', round(num_measurements_L1L1(i)));
    fprintf(fid2, '%d\n', round(num_measurements_Warn(i)));
    fprintf(fid3, '%3.3f\n', max(log10(reconstruction_error_L1L1(i)),-9)+9);
    fprintf(fid4, '%3.3f\n', max(log10(reconstruction_error_Warn(i)),-9)+9);
end
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); 
% =========================================================================


