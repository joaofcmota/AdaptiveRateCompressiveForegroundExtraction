function [xsol, output] = decoptL1L1Noise(b, sigma, w, beta, mode, arg1, arg2, ...
    varargin)

% [xsol, output] = decoptL1L1(b, sigma, w, beta, mode, arg1, arg2, varargin)
%
% Solves 
%
%             minimize     ||x||_1 + beta*||x-w||_1
%                x
%             subject to   ||A*x - b|| <= sigma
%
% where b: m x 1,  A: m x n,  w: n x 1, and beta > 0. There should hold 
% m > n. We use the algorithm in 
%
%  Quoc Tran-Dinh and Volkan Cevher, Optimal Primal-Dual Algorithms 
%  for Constrained Convex Minimization. Tech. Report, LIONS, EPFL (2014).
%
% This function has two modes: in mode 1, it uses an
% explicit matrix A; in mode 2, access to A, namely A*x and A'*y, is done 
% implicitly through function calls.
%
% If mode == 1, arg1 is the matrix A
%               arg2 is the pseudo inverse of A: pinv(A)
%
% If mode == 2, arg1 is a function handler to A*x
%               arg2 is a function handler to A'*y
%
% Inputs:
%   - b: m x 1 vector
%   - sigma: positive
%   - w: n x 1 vector
%   - beta: positive number. Theory in [1] recommends setting beta = 1
%   - mode: either 1 or 2
%   - arg1: in mode 1, an m x n matrix; in mode 2, a function handler
%   - arg2: in mode 2, an n x m matrix; in mode 2, a function handler 
%
% Optional input: maximum number of iterations (default: 5000)
%
% Outputs:
%   - x_opt: solution of (1)
%   - k: number of iterations

% =========================================================================


% =========================================================================
% Check inputs

m = length(b);
n = length(w);

if beta <= 0 || sigma <= 0
    error('beta and sigma should be positive');
end

if mode == 1
    A = arg1; 
elseif mode == 2
    A  = arg1;
    AT = arg2;
    AAT = @(x) A(AT(x));
else
    error('Mode not recognized');
end

if nargin > 8
    error('There should be only one optional argument')
end

if ~isempty(varargin)
    MAX_ITER = varargin{1};
else
    MAX_ITER = 5000;
end
% =========================================================================


% =========================================================================
% Parameters

tolx        = 1e-6;
verbosity   = 0;
% =========================================================================

% Set the parameters.
param = decoptParamSettings('MaxIters', MAX_ITER, 'RelTolX', tolx, ...
    'RelTolFeas', tolx, 'Algorithm', 3, 'Verbosity', verbosity);
param.PrintLength = 75;

% Check the inputs.
if nargin < 4, param = []; end
if nargin < 3, error('At least three inputs must be provided!'); end

preprocess_time = tic;
if mode == 1
    [pdata, opts] = decoptPreSolver('L1-L1NS', A, b, 'X0', w, 'NoiseLevel', sigma);
else
    [pdata, opts] = decoptPreSolver('L1-L1NS', A, b, 'NX', n, 'X0', w, ...
        'AT', AT, 'NoiseLevel', sigma);
end

pdata.w     = w;
pdata.beta  = beta;
%pdata.sigma = sigma;

% Print a message
% if param.Verbosity >= 0
%     fprintf([' + DECOPT is solving ', opts.pstrType, ' problem ...\n']);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameter for kicking the gap ...
resetKickingParameters();
kickPars.isOrthA = false;

% Evaluate the normAtA of the linear operator A.
normAtA = decoptNormAtAeval('operator', pdata.Aopr, param.PwMaxIters, ...
    param.PwRelTol, pdata.ATopr, pdata.nx);

% Define the proximal operators.
proxOpers = decoptProxForm(pdata, opts);

% Finish the preprocessing ...
preprocess_time = toc(preprocess_time);

% Call the internal solver.
solve_time = tic;
if param.Algorithm == 1
    
    alg_str        = 'Bregman smoothing 1P2D';
    [xsol, output] = decopt1P2DSolver(pdata.Aopr, pdata.ATopr, ...
        pdata.b, proxOpers, pdata.lbx, pdata.ubx, ...
        pdata.x0, normAtA, opts.weights, ...
        param, kickPars, pdata.varargin{:});
    
elseif param.Algorithm == 2
    
    alg_str        = 'Bregman smoothing 2P1D';
    [xsol, output] = decopt2P1DSolver(pdata.Aopr, pdata.ATopr, ...
        pdata.b, proxOpers, pdata.lbx, pdata.ubx, ...
        pdata.x0, normAtA, opts.weights, ...
        param, kickPars, pdata.varargin{:});
    
elseif param.Algorithm == 3
    
    alg_str        = 'Augmented Lagrangian smoothing 1P2D';
    [xsol, output] = decopt1P2DALSolver(pdata.Aopr, pdata.ATopr, ...
        pdata.b, proxOpers, pdata.lbx, pdata.ubx, ...
        pdata.x0, normAtA, opts.weights,  ...
        param, kickPars, pdata.varargin{:});
    
elseif param.Algorithm == 4
    
    alg_str        = 'Augmented Lagrangian smoothing 2P1D';
    [xsol, output] = decopt2P1DALSolver(pdata.Aopr, pdata.ATopr, ...
        pdata.b, proxOpers, pdata.lbx, pdata.ubx, ...
        pdata.x0, normAtA, opts.weights,  ...
        param, kickPars, pdata.varargin{:});
    
else
    error('This solver does not exist right now!');
end
solve_time = toc(solve_time);

% Evaluate the overall computational time.
output.alg        = alg_str;
output.pre_time   = preprocess_time;
output.solve_time = solve_time;
output.total_time = preprocess_time + solve_time;
output.param      = param;
%output.data       = pdata;

% Print output.
decoptPrintMessages(param, 'output', output, param.PrintLength);

% The nest function that resets the kicking parameters.
    function resetKickingParameters()
        
        % These parameters can be fixed for all the problems.
        kickPars.decrDel = 0.99;
        kickPars.lbGamma = 1e-5;     % Default: 1e-6;
        kickPars.maxGam  = 1.0e+100; % Default: 1e8.
        kickPars.minGam  = 1.0e-6;
        
        % These parameters can be changed for different problems.
        kickPars.gamFact = 0;        % default: 1.0
        kickPars.incrGam = 1.005;     % default: 1.05.
        kickPars.decrGam = 1.00;
        kickPars.gamma0  = [];
        
        % Change the value if they are provided.
        if isfield(opts, 'incrGam')
            if ~isempty(opts.incrGam), kickPars.incrGam = opts.incrGam; end
        end
    end

end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.
