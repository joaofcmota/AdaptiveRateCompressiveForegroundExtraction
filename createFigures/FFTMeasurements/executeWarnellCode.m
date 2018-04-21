function [] = executeWarnellCode(Execution_parameters, TOL_NONZERO)

% Wrapper to execute the code in
%
%    [1] Warnell et al. "Adaptive Rate Compressive Sensing for Background
%        Subtraction." http://arxiv.org/abs/1401.0583, 2014
%
% The core functions were taken from 
%
%    [2] http://garrettwarnell.com/ARCS-1.0.zip
%
% Input:
%
%  Execution_parameters:
%    FILENAME_DATA:    string with filename of dataset
%    FILENAME_RESULTS: string with desired filename for storing results
%    FILENAME_FRAMES:  string with desired filename for storing frames
%    nFrames:          number of frames that will be considered (first one
%                      should be the background image)
%    seed_num:         positive number to make replication of exp. possible
%    save_frames:      vector with frame numbers of which we will save the
%                      estimated foreground and the reconstructed frame
%
%  TOL_NONZERO: to remove noise from the image sequence by preprocessing it


% =========================================================================
% Parameters

DEBUG_OUTPUT = 1;   % If true, outputs some info

% Values taken from [2] for the sequences PETS, convoy, marker_cam 
sig_b    = 0.0157;     
tau      = 0.1; 

% s_buffer differs from sequence to sequence. In [2], it is
%   95  for PETS
%   123 for convoy
%   41  for marker_cam
s_buffer = 150;

SHOW_IMAGES = 0;       % If = 1, show images and plot results; 0, otherwise 
SAVE_DATA   = 1;       % If = 1, save data; 0, otherwise
% =========================================================================


% =========================================================================
% Read input and set parameters

try
    % Filename of the dataset
    FILENAME_DATA = Execution_parameters.FILENAME_DATA;
    
    % Filename to save the results
    FILENAME_RESULTS = Execution_parameters.FILENAME_RESULTS;
    
    % Filename to save the frames
    FILENAME_FRAMES = Execution_parameters.FILENAME_FRAMES;

    % Number of frames that will be used
    nFrames = Execution_parameters.nFrames;
      
    % Random seed number
    seed_num = Execution_parameters.seed_num;
    
    % Vector with frame numbers to be saved
    save_frames = Execution_parameters.save_frames;
    
catch err
    fprintf('Some of the parameters of the input of the function executeWarnellCode were not set correctly')
    return
end

% Set up paths
addpath('codeWarnell/src/CS/Decoder/')
addpath('codeWarnell/src/CS/Encoder/')
addpath('codeWarnell/src/CS/PhaseTransition/')


% Set the seed for generating random numbers
seed = RandStream('mcg16807','Seed',seed_num);
%RandStream.setGlobalStream(seed);
RandStream.setDefaultStream(seed);

vec = @(x) x(:);
% =========================================================================


% =========================================================================
% If saving data is active, create folder to data and another for images

image_folder_name   = 'results_img';
results_folder_name = 'results';
if SAVE_DATA 
    if ~isdir(results_folder_name)
        system(['mkdir ', results_folder_name]);
    end
    if ~isdir(image_folder_name)
        system(['mkdir ', image_folder_name]);
    end
end
% =========================================================================


% =========================================================================
% Read input: .mat file should have the following fields
%   nFrames: number of frames
%   siz:     1 x 2 array with the size of the images
%   imStack: siz(1) x siz(2) x nFrames matrix with raw images
%            imStack(:,:,1) should contain a clean background image 

load(FILENAME_DATA);            

nFrames = Execution_parameters.nFrames;   % To re-write the one in datafile


% **********************************************************
% Truncate differences in imStack
%
% Although this may be cheating, it is necessary to remove
% camera noise.

for fr = 2 : nFrames
    im_diff = imStack(:,:,fr) - imStack(:,:,1);
    im_diff_tr = im_diff;
    im_diff_tr(abs(im_diff) <= TOL_NONZERO) = 0;
    
    imStack(:,:,fr) = imStack(:,:,1) + im_diff_tr;   
    %figure(1);clf;imagesc(uint8(imStack(:,:,fr)));colormap(gray);drawnow
    %pause
end
imStack = im2double(uint8(imStack));
% **********************************************************

% =========================================================================


% =========================================================================
% Most of the following code is adapted from ARCSCV.m in [2]

% Maximum number of iterations to use in cs decoder
maxIt_cs = 1000;

%cross validation: accuracy, confidence, and number of measurements
eps_cv = 0.2;
xi_cv = 0.1;
r_cv = ceil(8*(1/eps_cv)^2*log(1/(2*xi_cv)));

% Parameters of the image
Im_hr = siz(1);
In_hr = siz(2);
N_hr  = Im_hr*In_hr;

% *************************************************************************
% Phase diagram and compressive sensing parameters

% Load Fourier phase transition
P = load('Fourier_10252013.mat');

rng_seed  = P.masterinfo.phi_seed;
bpdn_eps  = P.masterinfo.bpdn_eps;
spgl_opts = P.masterinfo.spgl_opts;

% Turn off spgl1 output
spgl_opts.verbosity = 0;

% Generate phase diagram
success_error = 1e-3;
pt_diagram = sum(P.pt_data<success_error,3)/size(P.pt_data,3);
% *************************************************************************

% *************************************************************************
% Estimate static component in compressive and cross-validation domains
if DEBUG_OUTPUT
    fprintf('Starting background estimation...\n');
end

% Define functions used to simulate \Phi and \Psi    
CS_fun = @(x)CS_FT(x,N_hr,rng_seed);
CV_fun = @(x)CS_Bernoulli(x,r_cv,rng_seed);

% Calculate full CS measurement vector for the high-resolution background
beta = CS_fun(vec(imStack(:,:,1)));
zeta = CV_fun(vec(imStack(:,:,1)));

if DEBUG_OUTPUT
    fprintf('Background estimation complete!\n\n');
end
% *************************************************************************

% *************************************************************************
% Test frame processing

if DEBUG_OUTPUT
    fprintf('Starting test frame processing...\n');
end

% Ground-truth quantities to record
s       = zeros(nFrames,1);
M_s     = zeros(nFrames,1);
sigma_s = zeros(nFrames,1);
sig_bs  = zeros(nFrames,1);

% Adaptive-rate quantities to record
hats          = zeros(nFrames+1,1);
M_hats        = zeros(nFrames,1);
sigma_hats    = zeros(nFrames,1);
errorbound_cv = zeros(nFrames,1);
error_cv      = zeros(nFrames,1);
kstar         = zeros(nFrames,1);

% Initialize sparsity estimate
hats(1) = s_buffer;
% *************************************************************************

% ==========================
% Variables to store results

% Number of measurements 
numMeasurements = zeros(nFrames,1);  

% Error between reconstructed image and real image (camera noise removed)
reconstruction_error = zeros(nFrames,1);
% ==========================

% *************************************************************************
% Loop through frames to process

for t = 1 : nFrames
    
    % Pre-sensing ground truth information
    f = vec(imStack(:,:,t)) - vec(imStack(:,:,1));
    
    s(t) = sum(abs(f) > tau);

    [f_sort, idx] = sort(abs(f),'descend');
    f_s                    = zeros(N_hr,1); 
    f_s(idx(1:s(t)))       = f(idx(1:s(t)));
    f_hats                 = zeros(N_hr,1); 
    f_hats(idx(1:hats(t))) = f(idx(1:hats(t)));

    sigma_s(t)    = norm(f - f_s);
    sigma_hats(t) = norm(f - f_hats);
    sig_bs(t)     = std(f(abs(f) < tau));
    
    % Determine number of compressive measurements to use for this frame
    % based on s, hats, and the phase diagram
    M_s(t) = PhaseTransitionLookup(s(t),N_hr, P.masterinfo.all_delta_pts, ...
        P.masterinfo.all_rho_pts, pt_diagram);
    
    M_s(t) = min(N_hr, max(1, M_s(t)));
    
    M_hats(t) = PhaseTransitionLookup(hats(t), N_hr, ...
        P.masterinfo.all_delta_pts,P.masterinfo.all_rho_pts,pt_diagram);
    
    M_hats(t) = min(N_hr, max(1, M_hats(t)));

    numMeasurements(t) = M_hats(t);   

    % Adaptive-rate compressive sensing, background subtraction, and
    % decoding
    
    CS_fun   = @(x)CS_FT(x, M_hats(t), rng_seed);
    CV_fun   = @(x)CS_Bernoulli(x, r_cv, rng_seed);
    CS_fun_T = @(x)CS_IFT(x, N_hr, rng_seed);
    
    SPG_CS_fun = @(x,mode) SPG_Interface(x, mode, CS_fun, CS_fun_T);
    
    y = CS_fun(vec(imStack(:,:,t)));
    
    beta_t = (sqrt(N_hr)/sqrt(M_hats(t)))*beta(1:M_hats(t));
        
    xi = y - beta_t;
    
    %spgl_opts.iterations = min(maxIt_cs, floor(M_hats(t)/5));
    spgl_opts.iterations = maxIt_cs;
    
    hatf = spg_bpdn(SPG_CS_fun, xi, bpdn_eps, spgl_opts);
        
    [hatf_sort, idx] = sort(abs(hatf),'descend');        
    hatf_hats = zeros(N_hr,1); 
    hatf_hats(idx(1:hats(t))) = hatf(idx(1:hats(t)));
    
    error_cv(t) = norm(f - hatf_hats);
    
    f_t  = 255*(imStack(:,:,1) + reshape(hatf_hats, siz(1), siz(2)));
    im_t = 255*(imStack(:,:,t));
    reconstruction_error_t = norm(f_t - im_t)/norm(im_t);
    reconstruction_error(t) = reconstruction_error_t;
    fprintf('Measurements used          = %d\n', M_hats(t));        
    fprintf('Reconstruction error       = %f  (frame %d)\n\n', reconstruction_error_t, t);        
    
    if sum(save_frames == t) > 0
        save([image_folder_name, '/', FILENAME_FRAMES, 'frame_', num2str(t), '.mat' ], ...
            'f_t');
    end    

    if SHOW_IMAGES
        figure(1);clf;
        subplot 121;
        imshow(uint8(im_t))
        title('Real image');
        subplot 122;
        imshow(uint8(f_t))
        title('Reconstructed image');
        drawnow
    end        

    % Use cross-validation procedure to predict next s using Gaussian pdf
    % approximation and hypothesis testing
    
    chi = CV_fun(vec(imStack(:,:,t)));
    gamma = chi - zeta;    
    errorbound_cv(t) = (1+eps_cv)*norm(CV_fun(hatf_hats) - gamma);
        
    mu0  = (N_hr - hats(t))*sig_b^2;
    var0 = 2*(N_hr-hats(t))*sig_b^4;

    k    = (hats(t)+1) : N_hr;
    muk  = (N_hr-k)*sig_b^2 + 1/3*( k-hats(t) )*( tau^2 + tau + 1 );
    
    vark = 1/9*( (k-hats(t)).^2 - (k-hats(t)) ) * ( tau^2 + tau+1 )^2 + ...
           1/5*( k-hats(t) )*( tau^4 + tau^3 + tau^2 + tau+1 ) + ...
           ( (N_hr-k).^2 + 2*(N_hr-k) )*sig_b^4 + ...
           2/3*( N_hr-k ).*( k-hats(t) )*( tau^2 + tau+1 )*sig_b^2 + ...
           - muk.^2;
    qval0  = pdf('normal', errorbound_cv(t)^2, mu0, sqrt(var0));
    qvalsk = pdf('normal', errorbound_cv(t)^2, muk, sqrt(vark));
    
    [qvalsk_max, qvalsk_max_idx] = max(qvalsk);
    
    if(qval0 >= qvalsk_max || errorbound_cv(t)^2 < mu0)

        kstar(t)  = 0;
        hats(t+1) = sum(abs(hatf_hats) > tau) + s_buffer;
        
    else        
        kstar(t)  = k(qvalsk_max_idx(1));
        hats(t+1) = kstar(t) + s_buffer;
    end
    
    % Make sure hats is within the interval [1 N_hr]
    hats(t+1) = max(1, min(hats(t+1), N_hr));
      
end

if DEBUG_OUTPUT
    fprintf('Test frame processing complete!\n\n');
end
% *************************************************************************

% =========================================================================


if SHOW_IMAGES    
    figure(2);clf;
    plot(numMeasurements,'b-');
    ylim([0, max(numMeasurements)]);
    xlabel('Frames')
    title('Number of measurements')
    drawnow
    
    figure(3);clf;
    semilogy(reconstruction_error, 'r-');
    xlabel('Frames')
    title('Relative error: ||xhat - x||/||x||');
    drawnow

end

if SAVE_DATA 
    save([results_folder_name, '/' FILENAME_RESULTS], ...
        'Execution_parameters', 'numMeasurements', 'reconstruction_error');
end

