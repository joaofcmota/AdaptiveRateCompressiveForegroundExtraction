% PURPOSE: This function tests the binary SVM problem.
%
% DATE: 21.04.2014.
%
% INPUTS:   file_name  : the name of the data file.
%           method     : proximal-gradient or proximal-Newton method.
%

%%
isExistTFOCS = 0;

% Input data.
rho_f     = 100.0;
%file_name = 'rcv1_test.binary';
%file_name = 'rcv1_train.binary';
%webspam_wc_normalized_unigram.svm';%'ijcnn1.t';%w2a.txt'; %colon-cancer';%a4a.txt'; 
%file_name = 'real-sim';
%file_name = 'news20.binary';
file_name = 'w4a.txt';
%file_name = 'a2a.t';
%file_name = 'a1a.t';
%file_name = 'ijcnn1.t';
%file_name = 'colon-cancer';

% If LIBSVM does not exist then use default mat data.
mat_file_name = 'w4a.mat';
fprintf('+ Read the data from file: %s ...\n', file_name);

if exist('svm_model_matlab.h')
    svm_path = which('svm_model_matlab.h');
    svm_path = svm_path(1:end-25);
    % Addpath to the libsvm toolbox.
    addpath([svm_path, 'matlab/']);
end
data_path = which('scopt_logistic_exam.m');
dirData   = [data_path(1:end-21), 'logistic_data/'];

% Read data from the data set.
if exist('libsvmread.c')
    [yhat, What] = libsvmread(fullfile(dirData, file_name));
    %save([dirData, mat_file_name], 'yhat', 'What');
else
    mdata = load([dirData, mat_file_name]);
    yhat  = mdata.yhat;
    What  = mdata.What;
end
if isempty(yhat) || isempty(What)
    error('Can not read the data from the file!');
end
bias = 0;

Wold = What;
bhat = -bias*yhat;

p = size(What,2);
N = size(What,1);
fprintf(' + We are solving problem: %s (p = %d, N = %d)\n', ...
         file_name, p, N);

% The penalty parameter.
rho   = rho_f/N; 

% Generate the initial point.
x0    = 0*ones(p, 1);

%% Call the decomposition algorithm.
tolx     = 1e-6;
maxiters = 3000;
isprint  = 2;

% Set the parameters.
param.MaxIters      = maxiters;
param.Verbosity     = isprint;
param.RelTolX       = tolx*1e-2;
param.RelTolFeas    = tolx*1e-2;
param.saveHistMode  = 0;
param.Algorithm     = 1;
param.InnerMaxIters = 50;
param.adaptStepSize = 0;

isL1Norm = true;
if isL1Norm, strpType = 'HingeL1'; else strpType = 'HingeL2'; end

% Call the solver.
[xsol_dc, output_dc] = decoptSolver(strpType, What, bias, param, ...
                                    'x0', x0, 'RegPar', rho, ...
                                    'Label', yhat);

% Evaluate the objective values.
Axb_dc  = What*xsol_dc - bhat;
fx_dc   = sum(max(1 - yhat.*Axb_dc, 0)) + rho*norm(xsol_dc, 1);
nnz_dc  = nnz(round(xsol_dc*1e10)*1e-10);

%% User-define prox-functions.
if isL1Norm
    proxOpers{1} = @(x, gamma, varargin)([proxL1norm(x(1:end-1), gamma); x(end)]);
    proxOpers{3} = @(x, varargin)( norm(rho.*x(1:end-1), 1) );
else
    proxOpers{1} = @(x, gamma, varargin)([proxSquareL2norm(x(1:end-1), gamma); x(end)]);
    proxOpers{3} = @(x, varargin)(0.5*sum(rho.*x(1:end-1).^2));
end
proxOpers{2} = @(x, gamma, varargin)(proxHingeLoss(x, gamma, yhat, yhat.^2));
proxOpers{4} = @(x, varargin)(sum(max(1 - yhat.*x, 0)));

% Regenerate the data.
What2 = [What, ones(N,1)];
bias2 = zeros(size(bias));
x02   = [x0; 0];
lbx   = [-inf*ones(p,1); 0];
ubx   = [+inf*ones(p,1); 10];

% Call the solver.
[xsol_dc2f, output_dc2] = decoptSolver('UserDef', What2, bias2, param, ...
                                       'x0', x02, 'RegPar', rho, ...
                                       'Label', yhat, 'Prox', proxOpers, ...
                                       'lbx', lbx, 'ubx', ubx);    
xsol_dc2 = xsol_dc2f(1:end-1);
new_bias = xsol_dc2f(end);

Axb_dc2  = What*xsol_dc2 - bhat;
fx_dc2   = sum(max(1 - yhat.*Axb_dc2, 0)) + rho*norm(xsol_dc2, 1);
nnz_dc2  = nnz(round(xsol_dc2*1e10)*1e-10);

% Print the outputs.
fprintf('**** FINAL OUTPUTS OF DECOMPOSITION ALGORITHM ****\n');
fprintf(' + Number of iterations: %5d\n', output_dc.iter);
fprintf(' + The objective value: %6.15f\n', fx_dc);
fprintf(' + CPU time: %6.4f\n', output_dc.total_time);
fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_dc);
fprintf('--------- USER-DEFINE-PROX-CASE -----------\n');
fprintf(' + Number of iterations: %5d\n', output_dc2.iter);
fprintf(' + The objective value: %6.15f\n', fx_dc2);
fprintf(' + CPU time: %6.4f\n', output_dc2.total_time);
fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_dc2);
fprintf('********************************************************\n');
                                
%%
%------------------------------------------------
% Compare with TFOCS. See http://cvxr.com/tfocs/.  
%------------------------------------------------
%%
if isExistTFOCS
    fprintf(' + CALL THE TFOCS SOLVER ...\n');
    Nx       = size(What, 1);
    nx       = size(What, 2);
    time1    = tic;
    affineF  = { diag(yhat)*What, bias*yhat; eye(nx), [] };
    if isL1Norm
        prox    = { prox_hingeDual(1, 1, -1), proj_linf(rho*Nx) };
    else
        prox    = { prox_hingeDual(1, 1, -1), proj_l2(rho*Nx) };
    end
    
    % Set options.
    opts.maxIts       = maxiters;
    opts.tol          = tolx;
    opts.alg          = 'N07';
    opts.restart      = -inf;

    % Call the TFOCS solver.
    [xsol_tfocs, output_tfocs, info] = tfocs_SCD([], affineF, prox, ...
                                                 1e-3, x0, [], opts);
    output_tfocs.time                = toc(time1);
    
    %% Evaluate the objective values.
    Axb_tfocs = What*xsol_tfocs - bhat;
    fx_tfocs  = sum(max(1 - yhat.*Axb_tfocs, 0)) + rho*norm(xsol_tfocs, 1);

    nnz_tfocs  = nnz(round(xsol_tfocs*1e10)*1e-10);
    fprintf('**** FINAL OUTPUTS OF TFOCS ****\n');
    fprintf(' + Number of iterations: %5d\n', output_tfocs.niter);
    fprintf(' + The objective value: %6.15f\n', fx_tfocs);
    fprintf(' + CPU time: %6.4f\n', output_tfocs.time);
    fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_tfocs);
    fprintf('********************************************************\n');
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF THE IMPLEMENTATION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%