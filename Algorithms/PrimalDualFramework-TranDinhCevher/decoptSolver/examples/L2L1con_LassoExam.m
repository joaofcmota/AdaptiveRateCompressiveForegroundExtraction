% TEST OUR SOLVER FOR THE L1-NORM CONSTRAINED LEAST SQUARES PROBLEM:
%                   
%                       min_x |A*x - b|_2  
%                       s.t. |x|_1 <= delta.
% Date: 31.12.2013
% Implemented by Quoc Tran-Dinh, LIONS, EPFL, Switzerland

%% Test other solvers.
isExistTFOCS = 1;
isPlotFigure = 1;

%% Problem size.
scale   = 1;
n       = round(scale*1000);
m       = round(scale*500); 
k       = round(scale*100); 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand('twister',0); randn('state',0);

% Noise level.
sigma   = 0.01;

% Generate matrix A & x.
cor_tau = 0.9;
if cor_tau > 0
    var0 = (1 - cor_tau)^2 / (1 - cor_tau^2); %initial variance
    A = zeros(m, n);
    A(:,1) = sqrt(var0)*randn(m, 1);
    for kk = 2:n
        A(:,kk) = cor_tau*A(:,kk-1) + (1 - cor_tau)*(randn(m,1));
    end
else
    A    = randn(m, n);
end

% Generate vector x_org.
x_org    = zeros(n, 1);
T        = randsample(n, k);
x_org(T) = randn(k, 1);

% Generate vector b.
b        = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0       = zeros(n, 1);

% Generate the constrained sparsity level.
delta    = 0.9*norm(x_org, 1);

%% Test the positive constrained L2/L1 problem.
tolx = 1e-6;

% Set the parameters.
param.MaxIters      = 3000;
param.Verbosity     = 2;
param.RelTolX       = tolx;
param.saveHistMode  = 0;
param.Algorithm     = 3;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;

% Call the solver.
[x1, out1] = decoptSolver('L2/L1con', A, b, param, 'NoiseLevel', ...
                          delta, 'x0', x0);

% Evaluate the objective values and feasibility gap.
fx1   = 0.5*norm(A*x1 - b, 2)^2; 
feas1 = max(norm(x1, 1) - delta, 0);

%% Plot the solution.
if isPlotFigure
	figure(1); title('The solutions');
    if isreal(x_org), stairs(x_org, 'g:*');  else stairs(abs(x_org), 'g:*'); end
    hold on;
    if isreal(x1), stairs(x1, 'r--o');  else stairs(abs(x1), 'r--o'); end
    shg;
end

%% Solve by TFOCS.
if isExistTFOCS
    % Optional parameters.
    opts          = [];
    opts.restart  = -inf;
    opts.alg      = 'N07';
    opts.tol      = param.RelTolX;
    opts.maxIts   = param.MaxIters;
    opts.countOps = 1;
    
    % Call the TFOCS solver.
    time2         = tic;
    [x2, out2, opts] = solver_LASSO(A, b, delta, [], opts);
    time2         = toc(time2);
    
    % Final solutions.
    out2.cntA  = out2.counts(end,3);
    out2.cntAt = out2.counts(end,2);
    out2.iter  = out2.niter;
    fx2        = 0.5*norm(A*x2 - b, 2)^2;
    feas2      = max( norm(x2, 1) - delta, 0);
end

%% Printing ...
fprintf('+ DECOM: L2/L1con - problem: f(x) = %3.7f, feas = %5.7f\n', fx1, feas1);
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistTFOCS
    fprintf('+ TFOCS: L2/L1con - problem: f(x) = %3.7f, feas = %5.7f\n', fx2, feas2);
    fprintf('+ TFOCS: Iterations: %4d, Time(s) = %3.4f\n', out2.iter, time2);
    fprintf('+ TFOCS: Number of Ax and ATy are %4d and %4d\n', out2.cntA, out2.cntAt);
    fprintf('+ TFOCS: Reconvery error: %4.7f\n', norm(x2 - x_org)/max(norm(x_org), 1));
end

%% Plot the figures.
if isPlotFigure && isExistTFOCS
    % Plot the solution.
    hold on; 
    if isreal(x2), stairs(x2, 'b-.s'); else stairs(abs(x2), 'b-.s'); end
    legend('Original signal', 'Decomp-solution', 'TFOCS-solution');
    hold off;
    shg;
end

%% END OF THE TEST.