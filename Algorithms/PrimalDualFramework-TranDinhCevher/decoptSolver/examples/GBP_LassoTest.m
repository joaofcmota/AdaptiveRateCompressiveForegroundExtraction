% TEST OUR SOLVER FOR THE BASIS PURSUIT PROBLEM OF THE FORM:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 31.12.2013

%% Test other solvers.
isExistYALL1  = 1;
isExistSGPL1  = 1; 

%plist         = [5:20, 5:20,21:25];
plist = [5, 6, 7, 8, 9, 10, 12, 15, 20, 25];
%plist = [5];

for pid = 1:length(plist)
    
%% Problem size.
scale   = plist(pid);
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
weights = ones(ng, 1);%rand(ng, 1) + 1;

% Noise level.
sigma   = 0*5e-2;

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
maxiters = 2000;
isprint  = 0;

% Set the parameters.
param.MaxIters      = maxiters;
param.Verbosity     = isprint;
param.RelTolX       = tolx;
param.RelTolFeas    = tolx;
param.saveHistMode  = 0;

%% Call the solver (Algorithm 1).
param.Algorithm     = 1;
param.adaptStepSize = 0;
time1a       = tic;
[x1a, out1a] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1a       = toc(time1a);

% Evaluate the objective values and feasibility gap.
fx1a   = sum(weights.*sqrt(mgroups*x1a.^2));
feas1a = norm(A*x1a - b, 2); 

%% Call the solver (Algorithm 2).
param.Algorithm     = 2;
param.adaptStepSize = 0;
time1b       = tic;
[x1b, out1b] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1b       = toc(time1b);

% Evaluate the objective values and feasibility gap.
fx1b   = sum(weights.*sqrt(mgroups*x1b.^2));
feas1b = norm(A*x1b - b, 2); 

%% Call the solver (Algorithm 3a).
param.Algorithm     = 3;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;
time1c       = tic;
[x1c, out1c] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1c       = toc(time1c);

% Evaluate the objective values and feasibility gap.
fx1c   = sum(weights.*sqrt(mgroups*x1c.^2));
feas1c = norm(A*x1c - b, 2); 

%% Call the solver (Algorithm 4a).
param.Algorithm     = 4;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;
time1d       = tic;
[x1d, out1d] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1d       = toc(time1d);

% Evaluate the objective values and feasibility gap.
fx1d   = sum(weights.*sqrt(mgroups*x1d.^2));
feas1d = norm(A*x1d - b, 2); 

%% Call the solver (Algorithm 3b).
param.Algorithm     = 3;
param.InnerMaxIters = 5;
param.adaptStepSize = 0;
time1e       = tic;
[x1e, out1e] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1e       = toc(time1e);

% Evaluate the objective values and feasibility gap.
fx1e   = sum(weights.*sqrt(mgroups*x1e.^2));
feas1e = norm(A*x1e - b, 2); 

%% Call the solver (Algorithm 4b).
param.Algorithm     = 4;
param.InnerMaxIters = 5;
param.adaptStepSize = 0;
time1f       = tic;
[x1f, out1f] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1f       = toc(time1f);

% Evaluate the objective values and feasibility gap.
fx1f   = sum(weights.*sqrt(mgroups*x1f.^2));
feas1f = norm(A*x1f - b, 2); 

%% If exist YALL1, then compare with it.
if isExistYALL1 
    fprintf('---> YALL1 ...\n');
    opts3.tol      = tolx;
    opts3.maxit    = maxiters;
    opts3.print    = isprint;
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
    opts4 = spgSetParms('verbosity',  isprint, ...
                        'weights',    weights, ...
                        'bpTol',      3*tolx, ...
                        'optTol',     tolx, ...
                        'decTol',     3*tolx, ...
                        'iterations', maxiters);
    sigma_b = sigma*norm(b);
    time4             = tic;
    [x4, R, G, info4] = spg_group(A, b, gindices, sigma_b, opts4);
    time4             = toc(time4);
    
    % Evaluate the objective values and feasibility gap.
    fx4   = sum(weights.*sqrt(mgroups*x4.^2));
    feas4 = norm(A*x4 - b, 2); 
end


%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM(1): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1a, feas1a/norm(b));
fprintf('+ DECOM(1): Iterations: %4d, Time(s) = %3.4f\n', out1a.iter, out1a.total_time);
fprintf('+ DECOM(1): Number of Ax and ATy are %4d and %4d\n', out1a.cntA, out1a.cntAt);
fprintf('+ DECOM(1): Reconvery error: %4.12f\n', norm(x1a - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(2): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1b, feas1b/norm(b));
fprintf('+ DECOM(2): Iterations: %4d, Time(s) = %3.4f\n', out1b.iter, out1b.total_time);
fprintf('+ DECOM(2): Number of Ax and ATy are %4d and %4d\n', out1b.cntA, out1b.cntAt);
fprintf('+ DECOM(2): Reconvery error: %4.12f\n', norm(x1b - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(3): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1c, feas1c/norm(b));
fprintf('+ DECOM(3): Iterations: %4d, Time(s) = %3.4f\n', out1c.iter, out1c.total_time);
fprintf('+ DECOM(3): Number of Ax and ATy are %4d and %4d\n', out1c.cntA, out1c.cntAt);
fprintf('+ DECOM(3): Reconvery error: %4.12f\n', norm(x1c - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(4): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1d, feas1d/norm(b));
fprintf('+ DECOM(4): Iterations: %4d, Time(s) = %3.4f\n', out1d.iter, out1d.total_time);
fprintf('+ DECOM(4): Number of Ax and ATy are %4d and %4d\n', out1d.cntA, out1d.cntAt);
fprintf('+ DECOM(4): Reconvery error: %4.12f\n', norm(x1d - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(5): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1e, feas1e/norm(b));
fprintf('+ DECOM(5): Iterations: %4d, Time(s) = %3.4f\n', out1e.iter, out1e.total_time);
fprintf('+ DECOM(5): Number of Ax and ATy are %4d and %4d\n', out1e.cntA, out1e.cntAt);
fprintf('+ DECOM(5): Reconvery error: %4.12f\n', norm(x1e - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(6): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1f, feas1f/norm(b));
fprintf('+ DECOM(6): Iterations: %4d, Time(s) = %3.4f\n', out1f.iter, out1f.total_time);
fprintf('+ DECOM(6): Number of Ax and ATy are %4d and %4d\n', out1f.cntA, out1f.cntAt);
fprintf('+ DECOM(6): Reconvery error: %4.12f\n', norm(x1f - x_org)/max(norm(x_org), 1));

if isExistYALL1 
    fprintf('+ YALL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx3, feas3/norm(b));
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
%    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.12f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end
if isExistSGPL1 
    fprintf('+ SPGL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx4, feas4/norm(b));
    fprintf('+ SPGL1: Iterations: %4d, Time(s) = %3.4f\n', info4.iter, info4.timeTotal);
    fprintf('+ SPGL1: Number of Ax and ATy are %4d and %4d\n', info4.nProdA, info4.nProdAt);
    fprintf('+ SPGL1: Reconvery error: %4.12f\n', norm(x4 - x_org)/max(norm(x_org), 1));
end

%% clear A;
save_file = ['p_5noiseless_', num2str(m), '_', num2str(n), '_', num2str(ng), '_', num2str(pid), '.mat'];
save(save_file); 

end

%% END OF THE TEST.