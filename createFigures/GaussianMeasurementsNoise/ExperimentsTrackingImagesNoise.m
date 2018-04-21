function [] = ExperimentsTrackingImagesNoise(Algorithm_parameters, Execution_parameters)

% Executes experiments in the paper...
%
% The input Algorithm_parameters and Execution_parameters are structures,
% containing the following fields:
%
% Algorithm_parameters:
%   SOLVER_L1L1:   either 'ADMM' or 'DECOPT'
%   MAX_ITER_L1L1: number of iterations of the L1-L1 solver
%   TOL_NONZERO:   tolerance for setting a nonzero pixel to 0 (e.g., 20)
%   alpha:         exponential moving average filter parameter (in [0,1])
%   delta:         oversampling factor (positive number)
%   epsilon:       oversampling factor (in [0,1]) for the noisy case
%   sigma:         noise parameter (positive number)
%   SI_BOOST:      we set w = SI_BOOST*w (reduces number of bad components)
%   MotionEstimation_Blocksize: size of each block in ME algorithm
%   MotionEstimation_Searchlim: another parameter of the ME algorithm
%   Initial_sparsity1: estimate of the sparsity of the 1st foreground*
%   Initial_sparsity2: estimate of the sparsity of the 2nd foreground*
%
% *************************************************************************
% *NOTE: Both Initial_sparsity1 and Initial_sparsity2 should be a positive
%        number. However, they can also be the string 'optimal'. In that
%        case, the sparsity of the respective foreground is computed from 
%        the actual image, i.e., no estimates are used in the 
%        initialization.
% *************************************************************************
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
    % Solver
    SOLVER_L1L1 = Algorithm_parameters.SOLVER_L1L1;
    
    % Maximum number of iterations of L1-L1 solver
    MAX_ITER_L1L1 = Algorithm_parameters.MAX_ITER_L1L1;
    
    % Tolerance for setting a nonzero pixel to 0
    TOL_NONZERO = Algorithm_parameters.TOL_NONZERO;
    
    % Exponential moving average filter parameter
    alpha = Algorithm_parameters.alpha;
    
    % Oversampling factor
    delta = Algorithm_parameters.delta;
    
    % Oversampling factor
    epsilon = Algorithm_parameters.epsilon;
    
    % Noise parameter
    sigma = Algorithm_parameters.sigma;
    
    % Side information boosting factor: we set w = SI_BOOST*w to reduce the
    % number of bad components; If no boosting, set SI_BOOST = 1;
    SI_BOOST = Algorithm_parameters.SI_BOOST;
    
    % Motion estimation parameters
    opts_ME = struct(...
        'BlockSize', {Algorithm_parameters.MotionEstimation_Blocksize}, ...
        'SearchLimit', {Algorithm_parameters.MotionEstimation_Searchlim}...
        );
    
    Initial_sparsity1 = Algorithm_parameters.Initial_sparsity1;
    Initial_sparsity2 = Algorithm_parameters.Initial_sparsity2;

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
    fprintf('Some of the parameters of the input of the function ExperimentsTrackingImages were not set properly')
    return
end

% Set up solvers paths
if strcmpi(SOLVER_L1L1, 'ADMM')
    addpath('../../../Algorithms/basisPursuitPlusL1/');
    solverL1L1_name = @basisPursuitPlusL1;

elseif strcmpi(SOLVER_L1L1, 'DECOPT')
    addpath('../../../Algorithms/PrimalDualFramework-TranDinhCevher/decoptSolver/');
    addpath('../../../Algorithms/PrimalDualFramework-TranDinhCevher/decoptSolver/functions/');
    addpath('../../../Algorithms/PrimalDualFramework-TranDinhCevher/decoptSolver/proxs/');
    solverL1L1_name = @decoptL1L1Noise;
else
    error('Unknown solver');
end


% Set the seed for generating random numbers
seed = RandStream('mcg16807','Seed',seed_num);
%RandStream.setGlobalStream(seed);
RandStream.setDefaultStream(seed);

vec = @(x) x(:);


SHOW_IMAGES = 0;       % If = 1, show images and plot results; 0, otherwise 
SAVE_DATA   = 1;       % If = 1, save data; 0, otherwise

USE_MYGAUSSIAN_OP = 0; % If 1, Gaussian matrices are generated with my 
                       % myGaussianOp
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


% **********************************************************
% Compute sparsity of foreground in frames 2 and 3; then,
% compute the CS bound 

% Stores sparsity of all images
sparsity_foreground_images = zeros(nFrames, 1);

diff_im1 = imStack(:,:,2) - imStack(:,:,1); 
diff_im2 = imStack(:,:,3) - imStack(:,:,1); 

sparsity_foreground_images(1) = sum(sum(diff_im1~=0));
sparsity_foreground_images(2) = sum(sum(diff_im2~=0));

if strcmp(Initial_sparsity1, 'optimal')
    s1 = sparsity_foreground_images(1);
else
    s1 = Initial_sparsity1;
end

if strcmp(Initial_sparsity2, 'optimal')
    s2 = sparsity_foreground_images(2);
else
    s2 = Initial_sparsity2;
end

CS_measurements_im1 = ceil((1/(1-epsilon)^2)*(2*s1*log(n/s1) + (7/5)*s1 + 1));
CS_measurements_im2 = ceil((1/(1-epsilon)^2)*(2*s2*log(n/s2) + (7/5)*s2 + 1));
% **********************************************************

% =========================================================================


%%
% =========================================================================
% Algorithm
%
% Reconstructs difference between first frame (background) and each frame

% *************************************************************************
% Initialization: for 2nd and 3rd frames, use standard CS; assume 
% background is known; assume sparsity is known

fprintf('CS measurements (first  frame) = %d\n', CS_measurements_im1);
fprintf('CS measurements (second frame) = %d\n', CS_measurements_im2);
fprintf('Initializing algorithm: solve 2 BPs...\n');

% % This matrix will be fixed. In subsequent frames, m_k measurements will
% % mean the first m_k rows of this matrix
% A = (1/sqrt(m_init))*randn(m_init,n);

% Background image
Im_backgr = imStack(:,:,1);

opts_spgl1 = spgSetParms('verbosity',0);   % Turn off verbose mode of spgl1

% =======
% Frame 2

if USE_MYGAUSSIAN_OP
    A = myGaussianOp(CS_measurements_im1, n);    
else
    Am = (1/sqrt(CS_measurements_im1))*randn(CS_measurements_im1, n);                        
    A = @(x, y) Am*x;
end
    

% Background measurements + noise
y_b = A(vec(Im_backgr), 1) + (sigma/sqrt(CS_measurements_im1))*randn(CS_measurements_im1,1);

% Take CS measurements of second frame
y_2 = A(vec(imStack(:,:,2)), 1);

% Reconstruct foreground
if USE_MYGAUSSIAN_OP    
    f_2_diff = spgl1(A, y_2 - y_b, 0 , sigma, vec(Im_backgr), opts_spgl1);
else
    f_2_diff = spgl1(Am, y_2 - y_b, 0 , sigma, vec(Im_backgr), opts_spgl1);
end

% Reconstruct image by adding foreground to the background
f_2 = Im_backgr + reshape(f_2_diff, siz(1), siz(2));
% =======

% =======
% Frame 3

if USE_MYGAUSSIAN_OP
    A = myGaussianOp(CS_measurements_im2, n);
else
    Am = (1/sqrt(CS_measurements_im2))*randn(CS_measurements_im2, n);
    A = @(x, y) Am*x;
end

% Background measurements + noise
y_b = A(vec(Im_backgr), 1) + (sigma/sqrt(CS_measurements_im2))*randn(CS_measurements_im2,1);

% Take CS measurements of third frame
y_3 = A(vec(imStack(:,:,3)), 1);

% Reconstruct foreground
if USE_MYGAUSSIAN_OP
    f_3_diff = spgl1(A, y_3 - y_b, 0 , sigma, vec(Im_backgr), opts_spgl1);
else
    f_3_diff = spgl1(Am, y_3 - y_b, 0 , sigma, vec(Im_backgr), opts_spgl1);
end

% Reconstruct image by adding foreground to the background
f_3 = Im_backgr + reshape(f_3_diff, siz(1), siz(2));
% =======
% *************************************************************************


% ==========================
% Variables to store results

% Number of L1L1 measurements if we knew the real signal
numMeasurementsL1L1_oracle = zeros(nFrames,1);  
numMeasurementsL1L1_oracle(1:2) = [CS_measurements_im1; CS_measurements_im2];

% Number of CS measurements if we knew the real signal
numMeasurementsCS_oracle = numMeasurementsL1L1_oracle;

% Error between side information and the true image difference
estimation_error = zeros(nFrames,1);

% Error between reconstructed image and real image (camera noise removed)
reconstruction_error = zeros(nFrames,1);

reconstruction_error(1) = norm(f_2 - imStack(:,:,2))/norm(imStack(:,:,2));
reconstruction_error(2) = norm(f_3 - imStack(:,:,3))/norm(imStack(:,:,3));
% ==========================


% Previous number of measurements
m_k_prev = CS_measurements_im2;

% Estimate on theoretical upper bound
L1L1_boundEstimate = numMeasurementsL1L1_oracle;   

% Number of measurements 
num_measurements = numMeasurementsL1L1_oracle;       

im_prev_prev = f_2;
im_prev      = f_3;

fprintf('Entering algorithm...\n');

%%
for kk = 3 : nFrames
    
    % *********************************************************************
    % Predict next frame using the last two reconstructed frames
    [MVx, MVy] = Motion_Estnew(im_prev_prev, im_prev, opts_ME);
    img_pred   = mceNew(im_prev, im_prev, MVx, MVy, opts_ME.SearchLimit);    
    img_pred = imresize(img_pred, siz);
    % *********************************************************************
        
    estimation_error(kk) = norm(img_pred - imStack(:,:,kk))...
        /norm(imStack(:,:,kk));
    
    % *********************************************************************
    % Update oracle quantites (i.e., knowing the real image)
    
    w = vec(img_pred - Im_backgr);            % side info
    
    x = vec(imStack(:,:,kk) - Im_backgr);     % real image difference
    
    % Truncate x and w
    w(abs(w) <= TOL_NONZERO) = 0;
    x(abs(x) <= TOL_NONZERO) = 0;
    
    true_foreground = imStack(:,:,kk) - Im_backgr;
    
    sparsity_foreground_images(kk) = sum(sum(true_foreground~=0));
    sk = sparsity_foreground_images(kk);
    
    w = SI_BOOST*w;
    
    h_bar = sum( (x>0).*(x>w) ) + sum( (x<0).*(x<w) );
    supp_x = (x~=0);
    xi    = sum(w(~supp_x) ~= x(~supp_x)) - sum(w(supp_x) == x(supp_x));
    
    m_theory = 2*h_bar*log(n/(sk + 0.5*xi)) + 7/5*(sk + 0.5*xi) + 1;
    m_theory = ceil(m_theory/(1-epsilon)^2);
    
    numMeasurementsL1L1_oracle(kk) = m_theory;
    
    numMeasurementsCS_oracle(kk) = (2*sk*log(n/sk) + (7/5)*sk + 1)/(1-epsilon)^2;

    fprintf('CS   measurements (oracle) = %d\n', numMeasurementsCS_oracle(kk));
    fprintf('L1L1 measurements (oracle) = %d\n', m_theory);
    
    if SHOW_IMAGES
        figure(1);clf;
        subplot 121;
        imshow(uint8(imStack(:,:,kk)))
        title('Real image');
        subplot 122;
        imshow(uint8(img_pred))
        title('Estimated image');
        drawnow
    end        
    
    % *********************************************************************
    
    % *********************************************************************
    % Decide the number of measurements
    
    % Estimate oracle number of measurements using exponential moving average
    L1L1_boundEstimate(kk) = (1-alpha)*L1L1_boundEstimate(kk-1) + alpha*m_k_prev;
    
    % Measurements we will take for this frame
    m_k = round(L1L1_boundEstimate(kk)*(1 + delta));
    
    num_measurements(kk) = m_k;
    
    fprintf('Measurements used          = %d \n', m_k);
    % *********************************************************************
        
    
    % Generate random matrix, take measurements of frame and of background
    if USE_MYGAUSSIAN_OP
        A_k = myGaussianOp(m_k, n);
    else
        Am = (1/sqrt(m_k))*randn(m_k,n);
        A_k = @(x, y) Am*x;
    end    
    y_k = A_k(vec(imStack(:,:,kk)), 1);
    y_b = A_k(vec(Im_backgr), 1) + (sigma/sqrt(m_k))*randn(m_k,1);
            
    % Reconstruct difference of frames using L1-L1 minimization
    if USE_MYGAUSSIAN_OP
        A_handler  = @(x) A_k(x, 1);
        AT_handler = @(y) A_k(y, 2);
        f_k_diff = solverL1L1_name(y_k - y_b, sigma, w, 1, 2, A_handler, AT_handler, ...
            MAX_ITER_L1L1);
    else
        A_k_pinv = pinv(Am);
        f_k_diff = solverL1L1_name(y_k - y_b, sigma, w, 1, 1, Am, A_k_pinv, ...
            MAX_ITER_L1L1);        
    end
    
    % Reconstruct image
    f_k = Im_backgr + reshape(f_k_diff, siz(1), siz(2));
    
    if sum(save_frames == kk) > 0
        % figure(100);clf;
        % imshow(uint8(img_pred));
        % drawnow;
        % saveas(gcf, [image_folder_name, '/', 'estimated_frame_', num2str(kk), '.pdf']);
        % figure(101);clf;
        % imshow(uint8(f_k));
        % drawnow;
        % saveas(gcf, [image_folder_name, '/', 'reconstructed_frame_', num2str(kk), '.pdf']);
        % figure(102);clf;
        % imshow(uint8(img_pred-imStack(:,:,1)));
        % drawnow;
        % saveas(gcf, [image_folder_name, '/', 'estimated_foreground_frame_', num2str(kk), '.pdf']);
        % figure(103);clf;
        % imshow(uint8(f_k-imStack(:,:,1)));
        % drawnow;
        % saveas(gcf, [image_folder_name, '/', 'reconstructed_foreground_frame_', num2str(kk), '.pdf']);
                
        save([image_folder_name, '/', FILENAME_FRAMES, 'frame_', num2str(kk), '.mat' ], ...
            'img_pred', 'f_k');
            
    end    
    
    % *********************************************************************
    % Update estimates of number of measurements and filter them
    
    w_prev = vec(img_pred - Im_backgr);
    x_prev = vec(f_k      - Im_backgr);
    
    % Truncate x_prev and w_prev
    w_prev(abs(w_prev) <= TOL_NONZERO) = 0;
    x_prev(abs(x_prev) <= TOL_NONZERO) = 0;
    
    w_prev = SI_BOOST*w_prev;
                                                
    h_bar_prev = sum( (x_prev>0).*(x_prev>w_prev) ) ...
        + sum( (x_prev<0).*(x_prev<w_prev) );

    supp_x_prev = (x_prev~=0);
    
    xi_prev = sum(w_prev(~supp_x_prev) ~= x_prev(~supp_x_prev)) ...
        - sum(w_prev(supp_x_prev) == x_prev(supp_x_prev));
    
    sk_prev = sum(sum(x_prev~=0));
       
    if h_bar_prev ~= 0                        
        m_k_prev = 2*h_bar_prev*log(n/(sk_prev + 0.5*xi_prev)) ...
            + 7/5*(sk_prev + 0.5*xi_prev) + 1;
        m_k_prev = ceil(m_k_prev/(1-epsilon)^2);
    else
        fprintf('No bad components in previous iteration: reducing number of measurements.\n');
        m_k_prev = round(m_k_prev/1.2);
    end        
    % *********************************************************************
    
    reconstruction_error_k = norm(f_k - imStack(:,:,kk))/norm(imStack(:,:,kk));
        
    reconstruction_error(kk) = reconstruction_error_k;
    
    fprintf('Reconstruction error       = %f  (frame %d)\n\n', reconstruction_error_k, kk);        
        
    im_prev_prev = im_prev;
    im_prev = f_k;    
end
% =========================================================================


%%

if SHOW_IMAGES
    
    figure(2);clf;
    plot(num_measurements,'bo-');
    hold on;
    plot(numMeasurementsL1L1_oracle,'ko-');
    plot(L1L1_boundEstimate,'go-');
    plot(numMeasurementsCS_oracle,'rs-');
    ylim([0,1.1*max([numMeasurementsCS_oracle; num_measurements])])
    xlabel('Frames')
    title('Number of measurements')
    legend('Measurements', 'L1L1 bound (oracle)', 'L1L1 bound estimate', 'CS measurements (oracle)')
    drawnow
    
    if SAVE_DATA 
        saveas(gcf, [image_folder_name, '/', FILENAME_RESULTS, 'meas.pdf']);
    end
    
    figure(3);clf;
    semilogy(estimation_error, 'ob-');
    hold on;
    semilogy(reconstruction_error, 'rs-');
    xlabel('Frames')
    title('Relative error: ||xhat - x||/||x||');
    legend('Estimation', 'Reconstruction')
    drawnow

    if SAVE_DATA 
        saveas(gcf, [image_folder_name, '/', FILENAME_RESULTS, 'error.pdf']);
    end
    
end

if SAVE_DATA 
    save([results_folder_name, '/' FILENAME_RESULTS], ...
        'Algorithm_parameters', 'Execution_parameters', 'n', ...
        'numMeasurementsL1L1_oracle', 'L1L1_boundEstimate', 'num_measurements', ...
        'numMeasurementsCS_oracle', 'estimation_error', 'reconstruction_error', ...
        'sparsity_foreground_images');
end


