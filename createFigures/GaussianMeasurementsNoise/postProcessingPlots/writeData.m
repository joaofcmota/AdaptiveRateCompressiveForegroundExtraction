
%%
% =========================================================================
% Experimental performance of the algorithms: file Experiments.m

% Hall
load('../results/ResultsHall128x128_allFrames.mat');
FilenameCSmeasurements_oracle ='HallCSmeasurements_oracle.dat';
FilenameL1L1meas_oracle       ='HallL1L1meas_oracle.dat';
FilenameL1L1meas_estim        ='HallL1L1meas_estim.dat';
FilenameMeasurements          ='HallMeasurements.dat';

% % PETS
% load('../results/ADMM/ResultsPETS_allFrames.mat');
% FilenameCSmeasurements_oracle ='PETSCSmeasurements_oracle.dat';
% FilenameL1L1meas_oracle       ='PETSL1L1meas_oracle.dat';
% FilenameL1L1meas_estim        ='PETSL1L1meas_estim.dat';
% FilenameMeasurements          ='PETSMeasurements.dat';

nFrames = Execution_parameters.nFrames;

fid1   = fopen(FilenameCSmeasurements_oracle,   'w');
fid2   = fopen(FilenameL1L1meas_oracle,         'w');
fid3   = fopen(FilenameL1L1meas_estim,          'w');
fid4   = fopen(FilenameMeasurements,            'w');

for i = 1 : nFrames %- 1 % For PETS, remove last image
    
    fprintf(fid1,   '%d %d\n', i, min(round(numMeasurementsCS_oracle(i)), 5900));  % For Hall
    %fprintf(fid1,   '%d %d\n', i, round(numMeasurementsCS_oracle(i)));
    fprintf(fid2,   '%d %d\n', i, round(numMeasurementsL1L1_oracle(i)));
    fprintf(fid3,   '%d %d\n', i, round(L1L1_boundEstimate(i)));
    fprintf(fid4,   '%d %d\n', i, round(num_measurements(i)));
    
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
% =========================================================================



%%
% =========================================================================
% Relative errors

% Hall
FilenameEstim ='Hall_EstimError.dat';
FilenameRec   ='Hall_RecError.dat';

% % PETS
% FilenameEstim ='PETSEstimError.dat';
% FilenameRec   ='PETSRecError.dat';

fid1   = fopen(FilenameEstim,   'w');
fid2   = fopen(FilenameRec,   'w');

for i = 1 : nFrames %- 1 % For PETS, remove last image

    fprintf(fid1, '%d %3.3f\n', i, max(log10(estimation_error(i)),-9)+9);
    fprintf(fid2, '%d %3.3f\n', i, max(log10(reconstruction_error(i)),-9)+9);
    
end

fclose(fid1);
fclose(fid2);