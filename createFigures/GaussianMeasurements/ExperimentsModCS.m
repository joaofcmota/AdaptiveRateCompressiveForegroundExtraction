function [] = ExperimentsModCS(Algorithm_parameters, Execution_parameters)

% Executes experiments for Mod-CS, which uses a fixed number of
% measurements.
%
% The input Algorithm_parameters and Execution_parameters are structs,
% containing the following fields:
%
% Algorithm_parameters:
%   MAX_ITER:      number of iterations of the Mod-CS solver
%   TOL_NONZERO:   tolerance for setting a nonzero pixel to 0 (e.g., 20)
%   NUM_MEAS:      number of measurements (positive integer smaller than n)
%   PERC_ENERGY:   number between 0 and 1; we estimate the support of x as
%                  the set of indices of the previous foreground frames
%                  that contain PERC_ENERGY of the coefficients' energy.
%
% Execution_parameters:
%   FILENAME_DATA:    string with filename of dataset
%   FILENAME_RESULTS: string with desired filename for storing results
%   FILENAME_FRAMES:  string with desired filename for storing frames
%   nFrames:          number of frames that will be considered (first one
%                     should be the background image)
%   seed_num:         positive number to make replication of exp. possible
%   save_frames:      vector with frame numbers of which we will save the
%                     estimated foreground and the reconstructed frame
                     

% =========================================================================
% Read input and set parameters

try
    % Maximum number of iterations of Mod-CS solver
    MAX_ITER = Algorithm_parameters.MAX_ITER;
    
    % Tolerance for setting a nonzero pixel to 0
    TOL_NONZERO = Algorithm_parameters.TOL_NONZERO;
    
    % Number of measurements
    NUM_MEAS = Algorithm_parameters.NUM_MEAS;
    
    % Number of measurements
    PERC_ENERGY = Algorithm_parameters.PERC_ENERGY;
        
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
    fprintf('Some of the parameters of the input of the function ExperimentsModCS were not set correctly')
    return
end

if PERC_ENERGY < 0 || PERC_ENERGY > 1
    error('PERC_ENERGY has to be a number between 0 and 1')
end

% Set up solver path
addpath('../../Algorithms/otherApproaches/modifiedCS-Vaswani/')

% Set the seed for generating random numbers
seed = RandStream('mcg16807','Seed',seed_num);
%RandStream.setGlobalStream(seed);
RandStream.setDefaultStream(seed);

vec = @(x) x(:);


SHOW_IMAGES = 0;       % If = 1, show images and plot results; 0, otherwise 
SAVE_DATA   = 1;       % If = 1, save data; 0, otherwise

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
% Read input: .mat file should have the following fields
%   nFrames: number of frames
%   siz:     1 x 2 array with the size of the images
%   imStack: siz(1) x siz(2) x nFrames matrix with raw images
%            imStack(:,:,1) should contain a clean background image 

load(FILENAME_DATA);            

nFrames = Execution_parameters.nFrames;   % To re-write the one in datafile

n = siz(1)*siz(2);

if mod(NUM_MEAS,1) || NUM_MEAS > n || NUM_MEAS <= 0
    error('NUM_MEAS should be a positive integer smaller than n')
end

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
% **********************************************************

% =========================================================================


%%
% =========================================================================
% Algorithm
%
% Reconstructs difference between first frame (background) and each frame

% *************************************************************************
% Initialization: for 1st frame, use standard CS; assume known background  

fprintf('Constructing measurement matrix and computing its pseudo-inverse...\n');
A = randn(NUM_MEAS, n);
A_pinv = pinv(A);

fprintf('\nInitializing: solve BP...\n');


% % This matrix will be fixed. In subsequent frames, m_k measurements will
% % mean the first m_k rows of this matrix
% A = (1/sqrt(m_init))*randn(m_init,n);

% Background image
Im_backgr = imStack(:,:,1);

opts_spgl1 = spgSetParms('verbosity',0);   % Turn off verbose mode of spgl1

% =======
% Frame 2

% Background measurements
y_b = A*vec(Im_backgr);

% Take CS measurements of second frame
y_2 = A*vec(imStack(:,:,2));

% Reconstruct foreground
f_2_diff = spgl1(A, y_2 - y_b, 0 , 0, vec(Im_backgr), opts_spgl1);

% Reconstruct image by adding foreground to the background
f_2 = Im_backgr + reshape(f_2_diff, siz(1), siz(2));
% =======

% ==========================
% Variables to store results

% Error between reconstructed image and real image (camera noise removed)
reconstruction_error_MOCS = zeros(nFrames,1);
reconstruction_error_MOCS(1) = norm(f_2 - imStack(:,:,2))/norm(imStack(:,:,2));
% ==========================

f_prev_diff = f_2_diff;

Difference_matrix_aux = ones(n);
Difference_matrix = tril(Difference_matrix_aux);

fprintf('Entering algorithm...\n');
%%
for kk = 2 : nFrames
    
    % Estimate the support of previous signal: where 90% of the 
    % coefficients' energy lies
    [amps, ind] = sort(abs(f_prev_diff), 'descend');
    total_energy = sum(amps.^2);
    cumulative_energy = Difference_matrix*(amps.^2);
    index_PERC_energy = find(cumulative_energy >= PERC_ENERGY*total_energy, 1);
    T = zeros(n,1);
    T(ind(1:index_PERC_energy)) = 1;
    T = (T~=0);
        
    % Take measurements of current frame
    y_k = A*vec(imStack(:,:,kk));
            
    % Reconstruct difference of frames using Modified-CS    
    [f_k_diff, iter_ADMM] = ModifiedCS(y_k - y_b, T, 1, A, A_pinv, MAX_ITER);         
    
    % Reconstruct image
    f_k = Im_backgr + reshape(f_k_diff, siz(1), siz(2));
       
    reconstruction_error_k = norm(f_k - imStack(:,:,kk))/norm(imStack(:,:,kk));
        
    reconstruction_error_MOCS(kk) = reconstruction_error_k;
    
    fprintf('Reconstruction error       = %f  (frame %d)\n\n', ...
        reconstruction_error_k, kk);
    f_prev_diff = f_k_diff;
end
% =========================================================================


%%

if SHOW_IMAGES
       
    figure(2);clf;        
    semilogy(reconstruction_error_MOCS, 'rs-');
    xlabel('Frames')
    title('Relative reconstruction error: ||xhat - x||/||x||');    
    drawnow

    if SAVE_DATA 
        saveas(gcf, [image_folder_name, '/', FILENAME_RESULTS, 'error.pdf']);
    end
    
end

if SAVE_DATA 
    save([results_folder_name, '/' FILENAME_RESULTS], ...
        'Algorithm_parameters', 'Execution_parameters', 'n', ...
        'reconstruction_error_MOCS');
end


