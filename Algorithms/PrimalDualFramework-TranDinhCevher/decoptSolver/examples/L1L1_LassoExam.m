% TEST OUR SOLVER FOR THE UNCONSTRAINED L1/L1 PROBLEM OF THE FORM:
%       
%                   min_x |A*x - b|_1 + |w.*x|_1
%
% Date: 31.12.2013
% Implemented by Quoc Tran-Dinh, LIONS, EPFL, Switzerland

%% Test other solvers.
isExistYALL1 = 1;
isPlotFigure = 1;

%% Problem size.
scale   = 1;
n       = scale*1000;
m       = scale*500; 
k       = scale*100; 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand('twister',0); 
%randn('state',0);

% Noise level.
sigma   = 1.0e-2;

% Generate matrix A.
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
    A       = randn(m, n);
end

% Generate vector b.
x_org    = zeros(n, 1);
T        = randsample(n, k);
x_org(T) = randn(k, 1);
noise    = thresh( randn(m, 1), 100);
b        = A*x_org + 2*sigma* noise;

% Generate an initial point.
x0       = zeros(n, 1);

% Generate the regularization parameter.
rho      = 0.5;

%% Test the unconstrained L1/L1 problem without weighted and no operators.
% Set the parameters.
param.MaxIters      = 5000;
param.Verbosity     = 2;   
param.RelTolX       = 1e-6;
param.saveHistMode  = 0;
param.Algorithm     = 3;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;


% Call the solver.
[x1, out1] = decoptSolver('L1/L1', A, b, param, 'RegPar', rho, 'x0', x0);
    
% Evaluate the objective values and feasibility gap.
fx1 = norm(A*x1 - b, 1) + norm(rho.*x1, 1);

%% If exist YALL1, then compare with it.
if isExistYALL1 
    opts3.tol      = param.RelTolX;
    opts3.maxit    = param.MaxIters;
    opts3.print    = param.Verbosity;
    opts3.nu       = rho;
    time3          = tic;
    [x3, out3]     = yall1(A, b, opts3);	
    time3          = toc(time3);
    % Evaluate the objective values and feasibility gap.
    fx3   = norm(A*x3 - b, 1) + norm(rho.*x3, 1);
end

%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: L1L1 - problem: f(x) = %3.7f\n', fx1);
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistYALL1 
    fprintf('+ YALL1: L1L1 - problem: f(x) = %3.7f\n', fx3);
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.7f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end

%% Plot the figures.
if isPlotFigure && isExistYALL1
    % Plot the solution.
    figure(1); title('The solutions');
    if isreal(x_org), stairs(x_org, 'g:*');  else stairs(abs(x_org), 'g:*'); end
    hold on;
    if isreal(x1), stairs(x1, 'r--o');  else stairs(abs(x1), 'r--o'); end
    hold on; 
    %if isreal(x2), stairs(x2, 'b-.s'); else stairs(abs(x2), 'b-.s'); end
    %legend('Original signal', 'Decomp-solution', 'TFOCS-solution');
    if isreal(x3), stairs(x3, 'b-.s'); else stairs(abs(x3), 'b-.s'); end
    legend('Original signal', 'Decomp-solution', 'YALL1');
    hold off;
    shg;
end

%% END OF THE TEST.