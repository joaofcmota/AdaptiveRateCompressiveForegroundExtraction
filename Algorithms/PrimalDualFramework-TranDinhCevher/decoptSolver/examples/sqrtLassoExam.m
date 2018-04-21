% TEST OUR SOLVER FOR THE UNCONSTRAINED SQUARE ROOT L1/L2 PROBLEM:
%
%                       min_x 0.5*|A*x - b|_2 + |w.*x|_1
% Date: 31.12.2013
% Implemented by Quoc Tran-Dinh, LIONS, EPFL, Switzerland

%% Test other solvers.
isExistTFOCS = 0;
isExistCvx   = 0;
isExistADMM  = 1;
isOperator   = 0;
isPlotFigure = 0;
isExistExactADMM = 1;

%% Problem size.
scale   = 10;
n       = round(scale*1000);
m       = round(scale*350); 
k       = round(scale*100); 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand('twister',0); randn('state',0);

% Noise level.
sigma   = 0.1;

% Generate matrix A.
cor_tau = 0.5;
if cor_tau > 0
  var0 = (1 - cor_tau)^2 / (1 - cor_tau^2); %initial variance
  A = zeros(m, n);
  A(:,1) = sqrt(var0)*randn(m, 1);
  for kk = 2:n
    A(:,kk) = cor_tau*A(:,kk-1) + (1 - cor_tau)*(randn(m,1));
  end
else
    A   = randn(m, n);
end

% Generate vector b.
x_org    = zeros(n, 1);
x_org(randperm(n, k)) = randn(k, 1);
b        = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0       = 0*ones(n, 1);

% Generate the regularization parameter.
c     = 1.1; 
ALPHA = 0.05; 
rho   = c*norminv(1 - ( ALPHA/(2*n) ), 0, sigma);
rho   = max(rho, 1e-4);

%% Test the unconstrained L1/L2 problem 

% Set the parameters.
param.MaxIters      = 5000;
param.Verbosity     = 2;    % Printing ...
param.RelTolX       = 1e-6;
param.RelTolFeas    = 1e-6;
param.saveHistMode  = 0;
param.Algorithm     = 1;
param.adaptStepSize = 0;
param.InnerMaxIters = 2;

% Call the solver.
[x1, out1] = decoptSolver('L1/sqrtL2', A, b, param, 'RegPar', rho, 'x0', x0);

% Evaluate the objective values and feasibility gap.
fx1 = norm(A*x1 - b, 2) + norm(rho.*x1, 1);

%% Plot the solution.
if 0%isPlotFigure
	figure(1); title('The solutions');
    if isreal(x_org), stairs(x_org, 'g:*');  else stairs(abs(x_org), 'g:*'); end
    hold on;
    if isreal(x1), stairs(x1, 'r--o');  else stairs(abs(x1), 'r--o'); end
    shg;
end

%% Test Inexact ADMM ...
if isExistADMM
    fprintf('+ We are running the inexact ADMM algorithm ...\n');
    
    % Set the options.
    opts = param;
    %opts.MaxIters = 1000;
    opts.PrintStep = 50;
    opts.RelTolX   = 1e-6;
    
    % Call the solver.
    time6      = tic;
    [x6, out6] = inexact_admm_sqrt_lasso(A, b, rho, opts);
    time6      = toc(time6);
    
    % Evaluate the objective values and feasibility gap.
    fx6 = norm(A*x6 - b, 2) + norm(rho.*x6, 1);
end

%% Test Exact ADMM ...
if isExistExactADMM
    fprintf('+ We are running the exact ADMM algorithm ...\n');
    
    % Set the options.
    opts = param;
    %opts.MaxIters = 1000;
    opts.PrintStep = 50;
    opts.RelTolX   = 1e-5;
    
    % Call the solver.
    time6b      = tic;
    [x6b, out6b] = exact_admm_sqrt_lasso(A, b, rho, opts);
    time6b      = toc(time6b);
    
    % Evaluate the objective values and feasibility gap.
    fx6b = norm(A*x6b - b, 2) + norm(rho.*x6b, 1);
end

%% Test with CVX.
if isExistCvx
    time4 = tic;
    if param.Verbosity <= 1, cvx_quiet true; end
    if param.Verbosity <= 0, cvx_quiet true; else cvx_quiet false; end
    cvx_begin;
        variable x_cvx(n, 1);
        minimize( norm(A*x_cvx - b, 2) + norm(rho.*x_cvx, 1) );
    cvx_end
    fx4   = cvx_optval;
    time4 = toc(time4);
end

%% If exist TFOCS, then compare with it.
if isExistTFOCS
    opts2           = tfocs_SCD; 
    %opts2.restart  = -inf;
    %opts2.alg      = 'N07';
    %opts2.alg      = 'GRA';
    %opts2.tol      = param.RelTolX;
    opts2.maxIts    = param.MaxIters;
    opts2.countOps  = 1;
    opts2.errFcn{1} = @(f, z, x)( norm(A*x - b, 2) + rho*norm(x, 1) );
    
    % Re define the regularize parameter.
    rho2           = c*sqrt(m)*norminv(1 - ( ALPHA/(2*n) ), 0, sigma);
    
    time2          = tic;
    %[x2, out2, info2] = tfocs_SCD(prox_l1(rho2/m), {A, -b}, ...
    %                              proj_l2(1/sqrt(m)), 1e-6, x0, [], opts2); 
    [x2, out2, info2] = tfocs_SCD(prox_l1(rho), {A, -b}, ...
                                  proj_l2(1.0), 1e-4, x0, [], opts2); 
    time2          = toc(time2);
    out2.cntA      = out2.counts(end,3);
    out2.cntAt     = out2.counts(end,2);
    out2.iter      = out2.niter;
    fx2            = norm(A*x2 - b, 2) + norm(rho.*x2, 1);
end

%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: sqrtLASSO-problem: f(x) = %3.7f\n', fx1);
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistCvx
    fprintf('+ CVX : sqrtLASSO-problem: f(x) = %3.7f\n', fx4);
    fprintf('+ CVX : Time(s) = %3.4f\n', time4);
    fprintf('+ CVX : Reconvery error: %4.7f\n', norm(x_cvx - x_org)/max(norm(x_org), 1));
end
if isExistTFOCS
    fprintf('+ TFOCS: sqrtLASSO-problem: f(x) = %3.7f\n', fx2);
    fprintf('+ TFOCS: Iterations: %4d, Time(s) = %3.4f\n', out2.iter, time2);
    fprintf('+ TFOCS: Number of Ax and ATy are %4d and %4d\n', out2.cntA, out2.cntAt);
    fprintf('+ TFOCS: Reconvery error: %4.7f\n', norm(x2 - x_org)/max(norm(x_org), 1));
end
if isExistADMM
    fprintf('+ iADMM: sqrtLASSO-problem: f(x) = %3.7f\n', fx6);
    fprintf('+ iADMM: Iterations: %4d, Time(s) = %3.4f\n', out6.iter, time6);
    fprintf('+ iADMM: Number of Ax and ATy are %4d and %4d\n', out6.cntA, out6.cntAt);
    fprintf('+ iADMM: Reconvery error: %4.7f\n', norm(x6 - x_org)/max(norm(x_org), 1));
end
if isExistExactADMM
    fprintf('+ eADMM: sqrtLASSO-problem: f(x) = %3.7f\n', fx6b);
    fprintf('+ eADMM: Iterations: %4d, Time(s) = %3.4f\n', out6b.iter, time6b);
    fprintf('+ eADMM: Number of Ax and ATy are %4d and %4d\n', out6b.cntA, out6b.cntAt);
    fprintf('+ eADMM: Reconvery error: %4.7f\n', norm(x6b - x_org)/max(norm(x_org), 1));
end

%% Plot the figures.
if 0%isPlotFigure && isExistADMM %&& isExistTFOCS
    % Plot the solution.
    hold on; 
    %if isreal(x2), stairs(x2, 'b-.s'); else stairs(abs(x2), 'b-.s'); end
    %legend('Original signal', 'Decomp-solution', 'TFOCS-solution');
    if isreal(x6), stairs(x6, 'b-.s'); else stairs(abs(x6), 'b-.s'); end
    legend('Original signal', 'Decomp-solution', 'iADMM-solution');
    hold off;
    shg;
end

%% Plot the objective values.
if isPlotFigure
    f1       = out1.hist.fx_val;
    f2       = out6.hist.fx_val;
    f2b      = out6b.hist.fx_val;
    best     = min([f1; f2; f2b]);
    if isExistCvx, best = min(best, fx4); end
    if isExistTFOCS, best = min(best, fx2); f3 = out2.err; end
    figure(randi([1,100], 1));
    semilogy((f1  - best)/abs(best), 'b');   hold on;
    semilogy((f2  - best)/abs(best), 'r');   hold on;
    semilogy((f2b - best)/abs(best), 'k--'); hold on;
    slg = {'DC', 'iADMM', 'eADMM'};
    if isExistTFOCS, semilogy((f3 - best)/abs(best), 'g--'); slg = {slg{:}, 'TFOCS'}; end
    legend(slg);
    shg;
end

%% END OF THE TEST.