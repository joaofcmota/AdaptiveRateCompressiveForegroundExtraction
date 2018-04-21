% PUROSE: The purpose of this code is to test some open-code solver for
% solving the following basis pursuit problem:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 16.06.2014.


addpath('../');
addpath('../functions');
addpath('../proxs');

%% Solvers for testing
isExistTFOCS  = 0;  % TFOCS can be downloaded at http://cvxr.com/tfocs/    
isExistYALL1  = 0;  % YALL1 can be downloaded at http://yall1.blogs.rice.edu
isExistSGPL1  = 1;  % SPGL1 can be downloaded at http://www.cs.ubc.ca/~mpf/spgl1/
isExistDecopt = 1;  % This is our code in the attachement.

%% Problem size.
scale   = 2; % Slcaling factor.
n       = scale*1000;
m       = scale*350; 
k       = scale*100; 

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
% Noise level.
sigma    = 1.0e-3;     % why are we using noise here? 
x_org    = zeros(n, 1);
T_aux    = randperm(n);   % ###
T = T_aux(1:k);           % ### 
x_org(T) = randn(k, 1);

% Generate matrix A & x.
A       = randn(m, n);

% Generate vector b.
b        = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0       = 0.*ones(n, 1);

%% Tolerance and maximum number of iterations.
tolx        = 1e-6;
maxiter     = 2000;
verbosity   = 2;

%% Test DECOPT
if isExistDecopt
    
    % Set the parameters.
    param = decoptParamSettings('MaxIters', maxiter, 'RelTolX', tolx, ...
            'RelTolFeas', tolx, 'Algorithm', 1, 'Verbosity', verbosity);
        
    % Call the solver.
    [x_decopt, out_decopt] = decoptSolver('BP', A, b, param, 'x0', x0);

    % Evaluate the objective values and feasibility gap.
    fx_decopt   = norm(x_decopt, 1);
    feas_decopt = norm(A*x_decopt - b, 2); 
end

%% If exist TFOCS, then compare with it.
if isExistTFOCS
    opts2             = [];
    opts2.restart     = -inf;
    opts2.tol         = tolx;
    opts2.maxIts      = maxiter;
    opts2.countOps    = 1;
    rho               = 1;
    mu                = .01*norm(x_org, Inf);
    opts2.errFcn{1}   = @(f, z, x)( norm(rho.*x, 1) );
    x02               = zeros(n, 1);
    time2             = tic;
    [x_tfocs, out_tfocs, info_tfocs] = solver_sBP(A, b, mu, x02, [], opts2);
    time2             = toc(time2);
    out_tfocs.cntA    = out_tfocs.counts(end,3);
    out_tfocs.cntAt   = out_tfocs.counts(end,2);
    out_tfocs.iter    = out_tfocs.niter;
    out_tfocs.time    = time2;
    
    % Evaluate the objective values and feasibility gap.
    fx_tfocs   = norm(x_tfocs, 1);
    feas_tfocs = norm(A*x_tfocs - b);
end

%% If exist YALL1, then compare with it.
if isExistYALL1 
    opts3.tol      = tolx;
    opts3.maxit    = maxiter;
    opts3.print    = verbosity;
    opts3.xs       = 1;
    time3          = tic;
    [x_yall1, out_yall1] = yall1(A, b, opts3);	
    time3          = toc(time3);
    out_yall1.time = time3;
    
    % Evaluate the objective values and feasibility gap.
    fx_yall1       = norm(x_yall1, 1);
    feas_yall1	   = norm(A*x_yall1 - b, 2); 
end

%% If exist YALL1, then compare with it.
if isExistSGPL1 
    opts4.bpTol      = tolx;
    opts4.optTol     = tolx;
    opts4.iterations = maxiter;
    opts4.maxMatvec  = 100000;
    opts4.verbosity  = verbosity;
    time4            = tic;
    [x_spgl1, R, G, info_spgl1]   = spg_bp(A, b, opts4);	
    
    time4            = toc(time4);
    out_spgl1        = info_spgl1;
    out_spgl1.time   = time4;
    
    % Evaluate the objective values and feasibility gap.
    fx_spgl1   = norm(x_spgl1, 1);
    feas_spgl1 = norm(A*x_spgl1 - b, 2); 
end


%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
if isExistDecopt
    fprintf('+ DECOM: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx_decopt, feas_decopt/norm(b));
    fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out_decopt.iter, out_decopt.total_time);
    fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out_decopt.cntA, out_decopt.cntAt);
    fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x_decopt - x_org)/max(norm(x_org), 1));
end
if isExistTFOCS
    fprintf('+ TFOCS: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx_tfocs, feas_tfocs/norm(b));
    fprintf('+ TFOCS: Iterations: %4d, Time(s) = %3.4f\n', out_tfocs.iter, out_tfocs.time);
    fprintf('+ TFOCS: Number of Ax and ATy are %4d and %4d\n', out_tfocs.cntA, out_tfocs.cntAt);
end
if isExistYALL1 
    fprintf('+ YALL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx_yall1, feas_yall1/norm(b));
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out_yall1.iter, out_yall1.time);
    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out_yall1.cntA, out_yall1.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.7f\n', norm(x_yall1 - x_org)/max(norm(x_org), 1));
end
if isExistSGPL1 
    fprintf('+ SPGL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx_spgl1, feas_spgl1/norm(b));
    fprintf('+ SPGL1: Iterations: %4d, Time(s) = %3.4f\n', out_spgl1.iter, out_spgl1.time);
    fprintf('+ SPGL1: Number of Ax and ATy are %4d and %4d\n', out_spgl1.nProdA, out_spgl1.nProdAt);
    fprintf('+ SPGL1: Reconvery error: %4.7f\n', norm(x_spgl1 - x_org)/max(norm(x_org), 1));
end

%% END OF THE TEST.