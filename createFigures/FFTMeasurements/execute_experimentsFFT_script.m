% Script to run the experiments

clear all
close all

% =========================================================================
% Add folders with data, and motion estimator to the path

addpath('../../datasets/')
addpath('SubME_1.6/');
% =========================================================================

% =========================================================================
% Hall sequence

fprintf('\n\n\nStart Hall sequence...\n\n');

TOL_NONZERO_HALL = 10;
NFRAMES_HALL     = 281;
SEED_NUM_HALL    = 1;
SAVE_FRAME_HALL  = [4,100,150,250];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_HALL}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'hall128x128.mat'}, ...
    'FILENAME_RESULTS', {'ResultsHall128x128_allFramesFFT.mat'}, ...
    'FILENAME_FRAMES', {'Hall_visualsFFT'}, ...
    'nFrames', {NFRAMES_HALL}, ...
    'seed_num', {SEED_NUM_HALL}, ...    
    'save_frames', {SAVE_FRAME_HALL} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'hall128x128.mat'}, ...
    'FILENAME_RESULTS', {'ResultsHall_allFramesFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'Hall_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_HALL}, ...
    'seed_num', {SEED_NUM_HALL}, ...
    'save_frames', {SAVE_FRAME_HALL} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_HALL);
% =========================================================================


% =========================================================================
% PETS sequence

fprintf('\n\n\nStart PETS sequence...\n\n');

TOL_NONZERO_PETS = 10;
NFRAMES_PETS     = 171;
SEED_NUM_PETS    = 1;
SAVE_FRAME_PETS  = [5,75,100,170];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_PETS}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'PETS09_S2L1.mat'}, ...
    'FILENAME_RESULTS', {'ResultsPETS_allFramesFFT.mat'}, ...
    'FILENAME_FRAMES', {'PETS_visualsFFT'}, ...
    'nFrames', {NFRAMES_PETS}, ...
    'seed_num', {SEED_NUM_PETS}, ...
    'save_frames', {SAVE_FRAME_PETS} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'PETS09_S2L1.mat'}, ...
    'FILENAME_RESULTS', {'ResultsPETS_allFramesFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'PETS_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_PETS}, ...
    'seed_num', {SEED_NUM_PETS}, ...
    'save_frames', {SAVE_FRAME_PETS} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_PETS);
% =========================================================================


% =========================================================================
% Stuttgart sequence

fprintf('\n\n\nStart Stuttgart sequence...\n\n');

TOL_NONZERO_STUT = 10;
NFRAMES_STUT     = 299;
SEED_NUM_STUT    = 1;
SAVE_FRAME_STUT  = [4,100,150,250];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_STUT}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'Stuttgart_144x192_0011_0310.mat'}, ...
    'FILENAME_RESULTS', {'Stuttgart_144x192_0011_0310FFT.mat'}, ...
    'FILENAME_FRAMES', {'Stuttgart_visualsFFT'}, ...
    'nFrames', {NFRAMES_STUT}, ...
    'seed_num', {SEED_NUM_STUT}, ...    
    'save_frames', {SAVE_FRAME_STUT} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'Stuttgart_144x192_0011_0310.mat'}, ...
    'FILENAME_RESULTS', {'Stuttgart_144x192_0011_0310FFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'Stuttgart_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_STUT}, ...
    'seed_num', {SEED_NUM_STUT}, ...
    'save_frames', {SAVE_FRAME_STUT} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_STUT);
% =========================================================================


% =========================================================================
% Vasco da Gama

fprintf('\n\n\nStart Vasco da Gama sequence...\n\n');

TOL_NONZERO_VCGM = 10;
NFRAMES_VCGM     = 891;
SEED_NUM_VCGM    = 1;
SAVE_FRAME_VCGM  = [4,200,600,850];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_VCGM}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'VascoDaGama_144x192_0344_1235.mat'}, ...
    'FILENAME_RESULTS', {'VascoDaGamaFFT.mat'}, ...
    'FILENAME_FRAMES', {'VascoDaGama_visualsFFT'}, ...
    'nFrames', {NFRAMES_VCGM}, ...
    'seed_num', {SEED_NUM_VCGM}, ...    
    'save_frames', {SAVE_FRAME_VCGM} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'VascoDaGama_144x192_0344_1235.mat'}, ...
    'FILENAME_RESULTS', {'VascoDaGamaFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'VascoDaGama_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_VCGM}, ...
    'seed_num', {SEED_NUM_VCGM}, ...
    'save_frames', {SAVE_FRAME_VCGM} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_VCGM);
% =========================================================================

%%
% =========================================================================
% baselineHighway sequence

fprintf('\n\n\nStart baselineHighway sequence...\n\n');

TOL_NONZERO_BSHG = 10;
NFRAMES_BSHG     = 499;
SEED_NUM_BSHG    = 1;
SAVE_FRAME_BSHG  = [10,100,250,450];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_BSHG}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'baselineHighway.mat'}, ...
    'FILENAME_RESULTS', {'ResultsBaselineHighwayFFT.mat'}, ...
    'FILENAME_FRAMES', {'BaselineHighwayFFT_visualsFFT'}, ...
    'nFrames', {NFRAMES_BSHG}, ...
    'seed_num', {SEED_NUM_BSHG}, ...    
    'save_frames', {SAVE_FRAME_BSHG} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'baselineHighway.mat'}, ...
    'FILENAME_RESULTS', {'ResultsBaselineHighwayFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'BaselineHighwayFFT_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_BSHG}, ...
    'seed_num', {SEED_NUM_BSHG}, ...
    'save_frames', {SAVE_FRAME_BSHG} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_BSHG);
% =========================================================================


% =========================================================================
% baselinePETS sequence

fprintf('\n\n\nStart baselinePETS sequence...\n\n');

TOL_NONZERO_BSPT = 10;
NFRAMES_BSPT     = 1199;
SEED_NUM_BSPT    = 1;
SAVE_FRAME_BSPT  = [100,500,750,1100];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_BSPT}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'baselinePETS2006.mat'}, ...
    'FILENAME_RESULTS', {'ResultsBaselinePetsFFT.mat'}, ...
    'FILENAME_FRAMES', {'BaselinePetsFFT_visualsFFT'}, ...
    'nFrames', {NFRAMES_BSPT}, ...
    'seed_num', {SEED_NUM_BSPT}, ...    
    'save_frames', {SAVE_FRAME_BSPT} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'baselinePETS2006.mat'}, ...
    'FILENAME_RESULTS', {'ResultsBaselinePetsFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'BaselinePetsFFT_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_BSPT}, ...
    'seed_num', {SEED_NUM_BSPT}, ...
    'save_frames', {SAVE_FRAME_BSPT} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_BSPT);
% =========================================================================


%%
% =========================================================================
% thermalPark sequence

fprintf('\n\n\nStart ThermalPark sequence...\n\n');

TOL_NONZERO_THPK = 10;
NFRAMES_THPK     = 30;%579;
SEED_NUM_THPK    = 1;
SAVE_FRAME_THPK  = [10,100,400,500];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_THPK}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'thermalPark.mat'}, ...
    'FILENAME_RESULTS', {'ResultsThermalParkFFT.mat'}, ...
    'FILENAME_FRAMES', {'ThermalParkFFT_visualsFFT'}, ...
    'nFrames', {NFRAMES_THPK}, ...
    'seed_num', {SEED_NUM_THPK}, ...    
    'save_frames', {SAVE_FRAME_THPK} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'thermalPark.mat'}, ...
    'FILENAME_RESULTS', {'ResultsThermalParkFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'ThermalParkFFT_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_THPK}, ...
    'seed_num', {SEED_NUM_THPK}, ...
    'save_frames', {SAVE_FRAME_THPK} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_THPK);
% =========================================================================


%%
% =========================================================================
% cameraJitterTraffic sequence

fprintf('\n\n\nStart cameraJitterTraffic sequence...\n\n');

TOL_NONZERO_CJTF = 10;
NFRAMES_CJTF     = 249;
SEED_NUM_CJTF    = 1;
SAVE_FRAME_CJTF  = [10,100,150,200];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_CJTF}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'cameraJitterTraffic.mat'}, ...
    'FILENAME_RESULTS', {'ResultsCameraJitterFFT.mat'}, ...
    'FILENAME_FRAMES', {'cameraJitterTrafficFFT_visualsFFT'}, ...
    'nFrames', {NFRAMES_CJTF}, ...
    'seed_num', {SEED_NUM_CJTF}, ...    
    'save_frames', {SAVE_FRAME_CJTF} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'cameraJitterTraffic.mat'}, ...
    'FILENAME_RESULTS', {'ResultsCameraJitterFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'cameraJitterTrafficFFT_Warnell'}, ...
    'nFrames', {NFRAMES_CJTF}, ...
    'seed_num', {SEED_NUM_CJTF}, ...
    'save_frames', {SAVE_FRAME_CJTF} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_CJTF);
% =========================================================================


% =========================================================================
% dynamicBackgroundCanoe sequence

fprintf('\n\n\nStart dynamicBackgroundCanoe sequence...\n\n');

TOL_NONZERO_CJTF = 10;
NFRAMES_CJTF     = 279;
SEED_NUM_CJTF    = 1;
SAVE_FRAME_CJTF  = [10,100,150,200];

Algorithm_parameters = struct('SOLVER_L1L1', {'ADMM'}, ...
    'MAX_ITER_L1L1', {1000}, ...
    'TOL_NONZERO', {TOL_NONZERO_CJTF}, ...
    'alpha', {0.5}, ...
    'delta', {0.01}, ...
    'SI_BOOST', {1.3}, ...
    'MotionEstimation_Blocksize', {8}, ...
    'MotionEstimation_Searchlim', {4}, ...
    'Initial_sparsity1', {150}, ...
    'Initial_sparsity2', {150} ...
    );

Execution_parameters = struct('FILENAME_DATA', {'dynamicBackgroundCanoe.mat'}, ...
    'FILENAME_RESULTS', {'ResultsDynamicBackgroundCanoeFFT.mat'}, ...
    'FILENAME_FRAMES', {'DynamicBackgroundCanoeFFT_visualsFFT'}, ...
    'nFrames', {NFRAMES_CJTF}, ...
    'seed_num', {SEED_NUM_CJTF}, ...    
    'save_frames', {SAVE_FRAME_CJTF} ...
    );

ExperimentsTrackingImagesFFT(Algorithm_parameters, Execution_parameters);

% Warnell CV method
Execution_parameters = struct('FILENAME_DATA', {'dynamicBackgroundCanoe.mat'}, ...
    'FILENAME_RESULTS', {'ResultsDynamicBackgroundCanoeFFT_Warnell.mat'}, ...
    'FILENAME_FRAMES', {'DynamicBackgroundCanoeFFT_visualsFFT_Warnell'}, ...
    'nFrames', {NFRAMES_CJTF}, ...
    'seed_num', {SEED_NUM_CJTF}, ...
    'save_frames', {SAVE_FRAME_CJTF} ...
    );

executeWarnellCode(Execution_parameters, TOL_NONZERO_CJTF);
% =========================================================================


fprintf('\n\n\nFinished!\n\n');
