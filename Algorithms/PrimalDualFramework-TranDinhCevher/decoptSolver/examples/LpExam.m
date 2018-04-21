% TEST OUR SOLVER FOR THE BASIS PURSUIT PROBLEM OF THE FORM:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 31.12.2013


%% Problem size.
K    = 10;
indn = randi([1, 50], K, 1);
N    = sum(indn);
M    = round(n/3);
s    = 0;
for i=1:K
    s1 = s + indn(i);
    index{i} = [s+1:s1];
    s = s1;
end

A    = randn(m, n);
c    = randn(n, 1);
xopt = rand(n, 1);
for pi=1:p
    xopti = xopt(index{pi});
    sp = sum(xopti);
    xopt(index{pi}) = xopti/sp;
end
b    = A*xopt;



% Print the problem size.
fprintf('+ The problem size [m, n, k] = [%d, %d, %d] ...\n', m, n, k);

%% Generate the input data.
%rand( 'twister', 0); randn('state',   0);

% Noise level.
sigma    = 1.0e-2;     % why are we using noise here? 
x_org    = zeros(n, 1);
T        = randsample(n, k);
x_org(T) = randn(k, 1);
x_org    = max(x_org, 0);
x_org    = x_org/sum(x_org);

% Generate matrix A and vector c.
A        = randn(m, n);
onev     = ones(1, n);
A        = [A; onev];
c        = randn(n, 1);

% Generate vector b.
b        = A*x_org;

% Generate an initial point.
x0       = 0.*ones(n, 1);

%% Test the BP problem.
tolx = 1e-6;
% Set the parameters.
param.MaxIters      = 3000;
param.Verbosity     = 2;
param.RelTolX       = tolx;
param.saveHistMode  = 0;
param.Algorithm     = 1;
param.InnerMaxIters = 2;
param.adaptStepSize = 0;

% Call the solver.
proxLpPos    = @(x, gamma, c)( max(0, x - gamma*c) );

% User-define proximal-functions.
proxOpers{1} = @(x, gamma, varargin)(proxLpPos(x, gamma, varargin{:}));
proxOpers{2} = @(x, gamma, varargin)(projL2norm(x, 1e-12));

proxOpers{3} = @(x, varargin)( varargin{1}'*x );
proxOpers{4} = @(x, varargin)(0);

% Call the solver with user-define prox-functions.
tic
[x1, out1] = decoptSolver('UserDef', A, b, param, 'x0', x0, 'Prox', proxOpers);
toc

% Evaluate the objective values and feasibility gap.


%% Call linprog to test.
isLinProg = 1;
if isLinProg
    tic
    opts       = optimset('Algorithm', 'Simplex');
    [x2 , fx2] = linprog(c, [], [], A, b, zeros(n,1), [], x0, opts);
    toc
    norm(x1 - x2)
end

%% END OF THE TEST.