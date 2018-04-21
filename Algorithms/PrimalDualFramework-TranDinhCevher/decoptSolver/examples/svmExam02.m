%%%%

%clear all;
%close all;
%clc

%train_list = {'a1a.txt', 'colon-cancer', 'rcv1_train.binary', 'news20.binary'};
%test_list  = {'a1a.txt', 'colon-cancer', 'rcv1_test.binary', 'news20.binary'};

train_list = { 'news20.binary' };
test_list = { '' };

% Data file.
%train_data_file = 'rcv1_train.binary'; 
%test_data_file  = 'rcv1_test.binary';

%train_data_file = 'news20.binary';
%test_data_file  = 'news20.binary';

% train_data_file = 'a2a.txt';
% test_data_file  = 'a2a.txt';

%train_data_file = 'colon-cancer';
%test_data_file  = 'colon-cancer';

for kkk = 1:1%length(train_list)
    
kkk2 = 1;    
train_data_file = train_list{kkk2};
test_data_file  = test_list{kkk2};

% If LIBSVM does not exist then use default mat data.
mat_file_name = 'w4a.mat';
fprintf('+ Read the data from file: %s ...\n', train_data_file);

if exist('svm_model_matlab.h')
    svm_path = which('svm_model_matlab.h');
    svm_path = svm_path(1:end-25);
    % Addpath to the libsvm toolbox.
    addpath([svm_path, 'matlab/']);
end
data_path = which('decoptSolver.m');
dirData   = [data_path(1:end-length('decoptSolver.m')), 'examples/svm/'];

% Read data from the data set.
if exist('libsvmread.c')
    [yhat_full, What_full] = libsvmread(fullfile(dirData, train_data_file));
    %save([dirData, mat_file_name], 'yhat', 'What');
else
    mdata = load([dirData, mat_file_name]);
    yhat_full  = mdata.yhat;
    What_full  = mdata.What;
end
if isempty(yhat_full) || isempty(What_full)
    error('Can not read the data from the file!');
end
bias = 0;

% If separate the data.
if isempty(test_data_file)
    full_size      = size(What_full, 1);
    train_size     = round(size(What_full, 1)*0.7);
    train_indices  = randperm(full_size, train_size);
    test_indices   = setdiff([1:full_size], train_indices);
    What           = What_full(train_indices, :);
    yhat           = yhat_full(train_indices, :);
    Wtest          = What_full(test_indices, :);
    ytest          = yhat_full(test_indices, :);
else
    What           = What_full;
    yhat           = yhat_full;
    [ytest, Wtest] = libsvmread(fullfile(dirData, test_data_file));
end

p = size(What,2);
N = size(What,1);
fprintf(' + We are solving problem: %s (p = %d, N = %d)\n', ...
          train_data_file, p, N);

% Generate the initial point.
x0 = 0*ones(p, 1);

%% Define the regularization parameters.

% Set of regularization parameters.
nPoints         = 10;
sC              = linspace(10^(-3), 10^3, nPoints);
number_of_Cs    = length(sC);
isLIBSVMRUN     = 1;

% Tracking along the parameter ...
for kk = 5:5%number_of_Cs

fprintf('%s\n', repmat('%', 1, 75));
fprintf('We are running the case: %d ...\n', kk);
    
% Take the regularization parameter.
C = sC(kk);

% Default parameters.
maxiters = 20000;
tolx     = 1e-13;

%%%%%%%%%%%%%%%%%%%%%%%%%%% THE DECOPT SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the parameters.
param.MaxIters      = maxiters;
param.Verbosity     = 0;
param.RelTolX       = tolx;
param.RelTolFeas    = tolx;
param.saveHistMode  = 0;
param.Algorithm     = 1;
param.InnerMaxIters = 50;
param.adaptStepSize = 1;

%% User-define prox-functions.
isL1Norm     = 0;
rho          = 1/C;
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
[N, p]  = size(What);
What2   = [What, ones(N,1)];
bias2   = zeros(size(bias));
x02     = [x0; 0];
lbx     = [-inf*ones(p,1); -1];
ubx     = [+inf*ones(p,1); +1];

%% Call the solver.
time_dc = tic;
%param.Verbosity = 2;
[xsol_dc2f, output_dc] = decoptSolver('UserDef', What2, bias2, param, ...
                                      'x0', x02, 'RegPar', rho, ...
                                      'Label', yhat, 'Prox', proxOpers, ...
                                      'GammaFactor', 1.006, ...
                                      'lbx', lbx, 'ubx', ubx);    
time_dc = toc(time_dc);

% Separate the solution ...                                  
xsol_dc2  = xsol_dc2f(1:end-1);
new_bias  = xsol_dc2f(end);
fsol.x    = xsol_dc2;
fsol.bias = new_bias;

% Evaluate the objective values.
[fx_dc, hx_dc, gx_dc, acc_dc]   = scopt_hinge_fx_eval(What, xsol_dc2, new_bias,... 
                                  C, yhat, isL1Norm);

[fx_dc2, hx_dc2, gx_dc2, acc_dc2] = scopt_hinge_fx_eval(Wtest, xsol_dc2, new_bias, ...
                                  C, ytest, isL1Norm);

nnz_dc  = nnz(round(xsol_dc2*1e10)*1e-10);
fprintf('**** FINAL OUTPUTS OF DECOMPOSITION ALGORITHM ****\n');
fprintf(' + Number of iterations: %5d\n', output_dc.iter);
fprintf(' + The objective value for test data: %6.15f\n', fx_dc);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*hx_dc, 0.5*gx_dc);
fprintf(' + Accuracy %1.5f\n',acc_dc);
fprintf(' + CPU time: %6.4f\n', output_dc.total_time);
fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_dc);
fprintf('---------- For the test data ----------\n');
fprintf(' + The objective value for test data: %6.15f\n', fx_dc2);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*hx_dc2, 0.5*gx_dc2);
fprintf(' + Accuracy %1.5f\n',acc_dc2);
fprintf('********************************************************\n');

fx_list1(:, kk)  = [fx_dc, fx_dc2];
gx_list1(:, kk)  = 0.5*[gx_dc, gx_dc2];
hx_list1(:, kk)  = C*[hx_dc, hx_dc2];
acc_list1(:, kk) = [acc_dc, acc_dc2];
time_list1(:,kk) = [time_dc];

%%%%%%%%%%%%%%%%%%%%%%%%%% CALL SVM SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Options for LibSVM ...
if isLIBSVMRUN
    fprintf('We are running the case: %d ...\n', kk);
    options_libsvm = ['-h 1 -s 0 -t 0 -e 0.001 -c ', num2str(C)];


    time_libsmv = tic;
    model       = svmtrain(yhat, What, options_libsvm);
    time_libsmv = toc(time_libsmv);

    %% Parameters from the model
    W_libsvm    = model.SVs' * model.sv_coef;
    bias_libsvm = -model.rho;

    [fx_svm, hx_svm, gx_svm, acc_svm] = scopt_hinge_fx_eval(What, ...
                                        W_libsvm,  bias_libsvm, C, yhat, false);

    [fx_svm2, hx_svm2, gx_svm2, acc_svm2] = scopt_hinge_fx_eval(Wtest, ...
                                            W_libsvm,  bias_libsvm, C, ytest, 0);

    fprintf(' ---- For the training data ----- ....\n');
    fprintf(' + The objective value for training data: %6.15f\n', fx_svm);
    fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*hx_svm, 0.5*gx_svm);
    fprintf(' + Accuracy %1.5f\n', acc_svm);
    fprintf(' + CPU time: %6.4f\n', time_libsmv);
    fprintf(' ---- For the test data ----- ....\n');
    fprintf(' + The objective value for test data: %6.15f\n', fx_svm2);
    fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*hx_svm2, 0.5*gx_svm2);
    fprintf(' + Accuracy %1.5f\n', acc_svm2);
end

    % Record the informations.
    if isLIBSVMRUN
        fx_list(:, kk)  = [fx_dc, fx_dc2, fx_svm, fx_svm2];
        gx_list(:, kk)  = 0.5*[gx_dc, gx_dc2, gx_svm, gx_svm2];
        hx_list(:, kk)  = C*[hx_dc, hx_dc2, hx_svm, hx_svm2];
        acc_list(:, kk) = [acc_dc, acc_dc2, acc_svm, acc_svm2];
        time_list(:,kk) = [time_dc, time_libsmv];
    end
end

%%
if ~isLIBSVMRUN
    fx_list(1:2,  :) = fx_list1;
	acc_list(1:2, :) = acc_list1;
    time_list(1, :)  = time_list1;
end

%% Plot if required.
isPlotFigure = 1;
if isPlotFigure
    figure('Name', 'Fx for training data');
    subplot(3, 2, 1);
    plot(fx_list(1,:), 'r--'); hold on;
    plot(fx_list(3,:), 'b:');  
    title('The objective values for the training data');
    legend('1P2D', 'LIBSVM');
    
    subplot(3, 2, 2);
    plot(fx_list(2,:), 'r--'); hold on;
    plot(fx_list(4,:), 'b:');  
    title('The objective values for the test data');
    legend('1P2D', 'LIBSVM');
    
    subplot(3, 2, 3);
    plot(acc_list(1,:), 'r--'); hold on;
    plot(acc_list(3,:), 'b:');  
    title('The accuracy values for the training data');
    legend('1P2D', 'LIBSVM');
    
    subplot(3, 2, 4);
    plot(acc_list(2,:), 'r--'); hold on;
    plot(acc_list(4,:), 'b:');  
    title('The accuracy values for the test data');
    legend('1P2D', 'LIBSVM');
    
    subplot(3, 2, 5);
    plot(time_list(1,:), 'r--'); hold on;
    plot(time_list(2,:), 'b:');  
    title('The cpu time');
    legend('1P2D', 'LIBSVM');
end

%%
save([dirData, train_data_file, '_new1_', num2str(kkk), '_', num2str(C), '.mat']);

end
%%%% END OF THE TEST ....