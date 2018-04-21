% TEST OUR SOLVER FOR THE BASIS PURSUIT PROBLEM OF THE FORM:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 31.12.2013

%% Test other solvers.
isExistTFOCS  = 0;
isExistYALL1  = 1;
isExistSGPL1  = 1; 
isPlotFigure  = 0;

%% Problem size.
scale   = 2;
n       = scale*1000;
m       = scale*350; 
k       = scale*100; 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand( 'twister', 0); randn('state',   0);

% Noise level.
sigma    = 1.0e-2;     % why are we using noise here? 
x_org    = zeros(n, 1);
T        = randsample(n, k);
x_org(T) = randn(k, 1);

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
b        = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0       = 0.*ones(n, 1);

%% Test the BP problem.
tolx = 1e-6;
% Set the parameters.
param.MaxIters      = 3000;
param.Verbosity     = 2;
param.RelTolX       = tolx;
param.saveHistMode  = 0;
param.Algorithm     = 5;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;

% Call the solver.
%[x1, out1] = decoptSolver('BP', A, b, param, 'x0', x0);

% User-define proximal-functions.
proxOpers{1} = @(x, gamma, varargin)(proxL1norm(x, gamma));
proxOpers{2} = @(x, gamma, varargin)(projL2norm(x, 1e-12));

proxOpers{3} = @(x, varargin)(norm(x(:), 1));
proxOpers{4} = @(x, varargin)(0);

% Call the solver with user-define prox-functions.
[x1, out1] = decoptSolver('UserDef', A, b, param, 'x0', x0, 'Prox', proxOpers, 'gamma0', 1e-2);

% Evaluate the objective values and feasibility gap.
fx1   = norm(x1, 1);
feas1 = norm(A*x1 - b, 2); 

%% If exist TFOCS, then compare with it.
if isExistTFOCS
    opts2             = [];
    opts2.restart     = -inf;
    opts2.tol         = tolx;
    opts2.maxIts      = param.MaxIters;
    opts2.countOps    = 1;
    opts2.errFcn{1}   = @(f, z, x)( norm(rho.*x, 1) );
    x02               = zeros(n, 1);
    time2             = tic;
    [x2, out2, info2] = solver_sBP(A, b, 1e-6, x02, [], opts2);
    time2             = toc(time2);
    out2.cntA         = out2.counts(end,3);
    out2.cntAt        = out2.counts(end,2);
    out2.iter         = out2.niter;
    
    % Evaluate the objective values and feasibility gap.
    fx2   = norm(x2, 1);
    feas2 = norm(A*x2 - b);
end

%% If exist YALL1, then compare with it.
if isExistYALL1 
    opts3.tol      = tolx;
    opts3.maxit    = param.MaxIters;
    opts3.print    = param.Verbosity;
    opts3.xs       = 1;
    time3          = tic;
    [x3, out3]     = yall1(A, b, opts3);	
    time3          = toc(time3);
    
    % Evaluate the objective values and feasibility gap.
    fx3            = norm(x3, 1);
    feas3          = norm(A*x3 - b, 2); 
end

%% If exist YALL1, then compare with it.
if isExistSGPL1 
    opts4.bpTol      = tolx;
    opts4.optTol     = tolx;
    opts4.iterations = param.MaxIters;
    opts4.maxMatvec  = 100000;
    opts4.verbosity  = 1;
    time4            = tic;
    [x4,R,G,info4]   = spg_bp(A, b, opts4);	
    
    time4         = toc(time4);
    
    % Evaluate the objective values and feasibility gap.
    fx4   = norm(x4, 1);
    feas4 = norm(A*x4 - b, 2); 
end


%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx1, feas1/norm(b));
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistTFOCS
    fprintf('+ TFOCS: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx2, feas2/norm(b));
    fprintf('+ TFOCS: Iterations: %4d, Time(s) = %3.4f\n', out2.iter, time2);
    fprintf('+ TFOCS: Number of Ax and ATy are %4d and %4d\n', out2.cntA, out2.cntAt);
end
if isExistYALL1 
    fprintf('+ YALL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx3, feas3/norm(b));
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.7f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end
if isExistSGPL1 
    fprintf('+ SPGL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx4, feas4/norm(b));
    fprintf('+ SPGL1: Iterations: %4d, Time(s) = %3.4f\n', info4.iter, info4.timeTotal);
    fprintf('+ SPGL1: Number of Ax and ATy are %4d and %4d\n', info4.nProdA, info4.nProdAt);
    fprintf('+ SPGL1: Reconvery error: %4.7f\n', norm(x4 - x_org)/max(norm(x_org), 1));
end

%% Plot the figures.
if 0 % isPlotFigure && isExistYALL1
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

%% Plot the objective values.
if isPlotFigure
    f1       = out1.hist.fx_val;
    best     = fx1;
    slg      = {'DC'};
    if isExistTFOCS, best = min(best, fx2); f2 = out2.err;     slg = {slg{:}, 'TFOCS'}; end
    if isExistYALL1, best = min(best, fx3); f3 = out3.objp;    slg = {slg{:}, 'YALL1'}; end
    if isExistSGPL1, best = min(best, fx4); f4 = info4.xNorm1; slg = {slg{:}, 'SPGL1'}; end
    figure(randi([1,100], 1));
    semilogy((f1 - best)/abs(best), 'b'); hold on;
    if isExistTFOCS, semilogy((f2 - best)/abs(best), 'r--'); hold on; end
    if isExistYALL1, semilogy((f3 - best)/abs(best), 'g.-'); hold on; end
    if isExistSGPL1, semilogy((f4 - best)/abs(best), 'k:'); hold  off; end
    legend(slg);
    shg;
end

%% END OF THE TEST.