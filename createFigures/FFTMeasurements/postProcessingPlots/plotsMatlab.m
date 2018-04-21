
%%
% =========================================================================
% Hall FFT

load('../results/ResultsHall128x128_allFramesFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsHall_allFramesFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,3500])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'HallFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'HallFFTError.png')
% =========================================================================


%%
% =========================================================================
% PETS FFT

load('../results/ResultsPETS_allFramesFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsPETS_allFramesFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

nFrames = Execution_parameters.nFrames;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,2500])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'PETSFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1(1:nFrames-1), 'b');
hold on;
semilogy(reconstruction_error_Warn(1:nFrames-1), 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'PETSFFTError.png')
% =========================================================================

%%
% =========================================================================
% Stuttgart FFT

load('../results/Stuttgart_144x192_0011_0310FFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/Stuttgart_144x192_0011_0310FFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,6000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'StutFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'StutFFTError.png')
% =========================================================================


%%
% =========================================================================
% VascoDaGamma FFT

load('../results/VascoDaGamaFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/VascoDaGamaFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,13000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'VGFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'VGFFTError.png')
% =========================================================================


%%
% =========================================================================
% baselineHighway FFT

load('../results/ResultsBaselineHighwayFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsBaselineHighwayFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,20000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'BSHGFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'BSHGFFTError.png')
% =========================================================================


%%
% =========================================================================
% baselinePETS FFT

load('../results/ResultsBaselinePetsFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsBaselinePetsFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,10000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'BSPTFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'BSPTFFTError.png')
% =========================================================================


%%
% =========================================================================
% thermalPark FFT

load('../results/ResultsThermalParkFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsThermalParkFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,4500])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'THPKFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'THPKFFTError.png')
% =========================================================================

%%
% =========================================================================
% cameraJitterTraffic FFT

load('../results/ResultsCameraJitterFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsCameraJitterFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,25000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'JTTRFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'JTTRFFTError.png')
% =========================================================================

%%
% =========================================================================
% dynamicBackgroundCanoe FFT

load('../results/ResultsDynamicBackgroundCanoeFFT.mat');

num_measurements_L1L1     = num_measurements;
reconstruction_error_L1L1 = reconstruction_error;

load('../results/ResultsDynamicBackgroundCanoeFFT_Warnell.mat');

num_measurements_Warn     = numMeasurements;
reconstruction_error_Warn = reconstruction_error;

figure(1);clf;
plot(numMeasurementsL1L1_oracle, 'g-.');
hold on;
plot(num_measurements_L1L1, 'b');
plot(num_measurements_Warn, 'r');
ylim([0,25000])
title('Number of measurements')
xlabel('Frame')
legend('L1-L1 oracle', 'L1-L1', 'Warnell')
saveas(gcf, 'DBCNFFTMeas.png')

figure(2);clf;
semilogy(reconstruction_error_L1L1, 'b');
hold on;
semilogy(reconstruction_error_Warn, 'r');
title('Relative error')
xlabel('Frame')
legend('L1-L1', 'Warnell')
saveas(gcf, 'DBCNFFTError.png')
% =========================================================================
