% TEST OUR SOLVER FOR THE BASIS PURSUIT PROBLEM OF THE FORM:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 31.12.2013

%% Test other solvers.
isExistYALL1  = 1;
isExistSGPL1  = 0; 

%% Problem size.
scale   = 5;
n       = scale*1024;
m       = scale*341; 
ng      = round(m/4);
K       = round(ng/8);

% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, K);

%% Generate the input data.
%rand( 'twister', 0); randn('state',   0);

% Generate the the group.
gindices = [];
while (length(unique(gindices)) ~= ng)
    gindices  = ceil(rand(n,1) * ng);
end
groups.ng = ng;
groups.indices = gindices;
mgroups = sparse(ng, n);
for i=1:ng
    mgroups(i, gindices == i) = 1;
end

% Generate the weight of the groups.
weights = rand(ng, 1) + 1;

% Noise level.
sigma      = 0*1e-3;

% Create K-group-sparse vector x0 and observation vector b
p          = randperm(ng); 
p          = p(1:K);
idx        = ismember(gindices, p);
nNz        = sum(idx);  % number of nonzeros
x_org      = zeros(n, 1); 
x_org(idx) = randn(sum(idx), 1);


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
    A    = randn(m, n);
end

% Generate vector b.
b        = A*x_org + sigma*randn(m, 1);

% Generate an initial point.
x0       = 0.*ones(n, 1);

%% Test the BP problem.
tolx     = 1e-6;
maxiters = 1000;

% Set the parameters.
param.MaxIters      = maxiters;
param.Verbosity     = 2;
param.RelTolX       = tolx;
param.RelTolFeas    = tolx;
param.saveHistMode  = 4;
param.Algorithm     = 4;
param.InnerMaxIters = 10;
param.InnerRelTol   = tolx;
param.adaptStepSize = 0;

% Call the solver.
[x1, out1] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups, 'RegPar', weights);

% Evaluate the objective values and feasibility gap.
fx1   = sum(weights.*sqrt(mgroups*x1.^2));
feas1 = norm(A*x1 - b, 2); 

%% If exist YALL1, then compare with it.
if isExistYALL1 
    fprintf('---> YALL1 ...\n');
    opts3.tol      = tolx;
    opts3.maxit    = maxiters;
    opts3.print    = param.Verbosity;
    opts3.xs       = 1;
    time3          = tic;
    
    %% Group YALL1 Solver.
    [x3, out3] = YALL1_group(A, b, groups.indices,...
                                   'StopTolerance', opts3.tol, ...
                                   'GrpWeights', weights, ...
                                   'overlap', false, ...
                                   'nonorthA', 1, ...
                                   'ExactLinSolve', 1, ...
                                   'Solver', 1, ...             % 1 - primal; 2 - dual
                                   'maxIter', opts3.maxit);
    time3 = toc(time3);
 
    % Evaluate the objective values and feasibility gap.
    fx3   = sum(weights.*sqrt(mgroups*x3.^2));
    feas3 = norm(A*x3 - b, 2); 
end

%% If exist YALL1, then compare with it.
if isExistSGPL1 
    fprintf('---> SPGL1 ...\n');
    opts4 = spgSetParms('verbosity', 1, ...
                        'weights',   weights, ...
                        'bpTol',     tolx, ...
                        'optTol',    5*tolx, ...
                        'decTol',    tolx, ...
                        'iterations',maxiters);
    sigma_b = sigma*norm(b);
    time4             = tic;
    [x4, R, G, info4] = spg_group(A, b, gindices, sigma_b, opts4);
    time4         = toc(time4);
    
    % Evaluate the objective values and feasibility gap.
    fx4   = sum(weights.*sqrt(mgroups*x4.^2));
    feas4 = norm(A*x4 - b, 2); 
end


%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1, feas1/norm(b));
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));
if isExistYALL1 
    fprintf('+ YALL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx3, feas3/norm(b));
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
%    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.7f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end
if isExistSGPL1 
    fprintf('+ SPGL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx4, feas4/norm(b));
    fprintf('+ SPGL1: Iterations: %4d, Time(s) = %3.4f\n', info4.iter, info4.timeTotal);
    fprintf('+ SPGL1: Number of Ax and ATy are %4d and %4d\n', info4.nProdA, info4.nProdAt);
    fprintf('+ SPGL1: Reconvery error: %4.7f\n', norm(x4 - x_org)/max(norm(x_org), 1));
end

%% END OF THE TEST.