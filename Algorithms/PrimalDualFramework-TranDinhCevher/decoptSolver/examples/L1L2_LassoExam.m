% TEST OUR SOLVER FOR THE UNCONSTRAINED L1/L2 PROBLEM OF THE FORM:
%
%                    min_x 0.5*|A*x - b|_2^2 + |w.*x|_1
% Date: 31.12.2013
% Create by Quoc Tran-Dinh, LIONS, EPFL, Switzerland.

%% Test other solvers.
isExistTFOCS = 0;
isExistYALL1 = 1;
isPlotFigure = 0;

%% Problem size.
scale   = 1;
n       = round(scale*1000);
m       = round(scale*400); 
k       = round(scale*100); 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand('twister',0); randn('state',0);

% Noise level.
sigma   = 0.01;

% Generate matrix A & x.
cor_tau = 0.5;
if cor_tau > 0
    var0 = (1 - cor_tau)^2 / (1 - cor_tau^2); %initial variance
    A = zeros(m, n);
    A(:,1) = sqrt(var0)*randn(m, 1);
    for kk = 2:n
        A(:,kk) = cor_tau*A(:,kk-1) + (1 - cor_tau)*(randn(m,1));
    end
else
    A = randn(m, n);
end

% Generate vector x_org.
x_org    = zeros(n, 1);
T        = randsample(n, k);
x_org(T) = randn(k, 1);

% Generate vector b.
b  = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0 = 0.*ones(n, 1);

% The regularization parameter.
rho = 0.1;    

%% Test the BP problem.
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
[x1, out1] = decoptSolver('L1/L2', A, b, param, 'RegPar', rho, 'x0', x0);

% Evaluate the objective values and feasibility gap.
fx1 = 0.5*norm(A*x1 - b, 2)^2 + norm(rho.*x1, 1);

%% If exist TFOCS, then compare with it.
if isExistTFOCS
    opts2          = [];
    opts2.restart  = -inf;
    opts2.alg      = 'N07';
    opts2.tol      = param.RelTolX;
    opts2.maxIts   = param.MaxIters;
    opts2.countOps = 1;
    x02            = zeros(n, 1);
    time2          = tic;
    [x2, out2, info2] = solver_L1RLS(A, b, rho, x02, opts2);
    time2          = toc(time2);
    out2.cntA  = out2.counts(end,3);
    out2.cntAt = out2.counts(end,2);
    out2.iter  = out2.niter;
    fx2        = 0.5*norm(A*x2 - b, 2)^2 + norm(rho.*x2, 1);
end

%% If exist YALL1, then compare with it.
if isExistYALL1 
    opts3          = [];
    opts3.tol      = param.RelTolX;
    opts3.maxit    = param.MaxIters;
    opts3.print    = param.Verbosity;
    opts3.rho      = 1;
    opts3.weights  = rho;
    opts3.xs       = 1;
    time3          = tic;
    [x3, out3]     = yall1(A, b, opts3);	
    time3          = toc(time3);
    % Evaluate the objective values and feasibility gap.
    fx3   = 0.5*norm(A*x3 - b, 2)^2 + norm(rho.*x3, 1);
end

%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: L1L2-problem: f(x) = %3.12f\n', fx1);
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistTFOCS
    fprintf('+ TFOCS: L1L2-problem: f(x) = %3.12f\n', fx2);
    fprintf('+ TFOCS: Iterations: %4d, Time(s) = %3.4f\n', out2.iter, time2);
    fprintf('+ TFOCS: Number of Ax and ATy are %4d and %4d\n', out2.cntA, out2.cntAt);
    fprintf('+ TFOCS: Reconvery error: %4.7f\n', norm(x2 - x_org)/max(norm(x_org), 1));
end
if isExistYALL1 
    fprintf('+ YALL1: L1L2-problem: f(x) = %3.12f\n', fx3);
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.7f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end

%% Plot the objective values.
if isPlotFigure
    f1       = out1.hist.fx_val;
    best     = fx1;
    if isExistYALL1, best = min(best, fx3); f2 = out3.objp(:); end
    figure(randi([1,100], 1));
    semilogy((f1 - best)/abs(best), 'b'); hold on;
    if isExistYALL1
        semilogy((f2 - best)/abs(best), 'r'); hold off;
    end
    legend('DC', 'YALL1');
    shg;
end
%% END OF THE TEST.