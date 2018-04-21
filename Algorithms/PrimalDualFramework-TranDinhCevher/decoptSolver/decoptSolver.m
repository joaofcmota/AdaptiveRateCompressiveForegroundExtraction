% FUNCTION: [xsol, output] = decoptSolver(pType, Aoper, b, param, varargin)
%
% PURPOSE: An implementation of the new decomposition algorithm for solving
%       one of the following problems:
%
%       1) L1/L2: The unconstrained L1/L2 LASSO problem:   
%                   min_x 0.5*|A*x - b|_2^2 + |W.*x|_1 
%
%       2) L1/L1: The L1/L1 problem:
%                       min_x |A*x - b|_1 + |W.*x|_1 
%
%       3) Square-root L1/L2 LASSO problem of the form:
%               min_x |A*x - b|_2 + |W.*x|
%
%       4) L1/L2con: The basis pursuit denoise (BPDN) problem:
%               min_x |W.*x|_1  s.t. |A*x - b|_2 <= delta
%
%       5) L2/L1con: The L1-constrained LASSO problem:
%               min_x 0.5*|A*x - b|_2^2 s.t. |x|_1 <= delta,
%
%       6) BP: The basis pursuit (BP) problem:
%                       min_x |W.*x|_1 s.t. A*x == b.
%
%       where A in C^{m x n} and b in C^m, m << n;
%             W is an n-positive regularization parameter vector
%               (it can be a positive scalar);
%             delta > 0 is a level noise;
%             lbx is the lower bound vector (if provided);    
%             ubx is the upper bound vector (if provided).
%
% REFERENCES: 
%  [1]. Quoc Tran-Dinh and Volkan Cevher, Optimal Primal-Dual Algorithms 
%       for Constrained Convex Minimization. Tech. Report, LIONS, EPFL (2014).
%
% INPUTS:
%   pType      - a string indicating type of problem. It can be one of the 
%                following sign: 
%                   'L1/L2', 'L1/L1', 'L1/sqrtL2', 'L1/L2con', 'L2/L1con' 
%                   and 'BP'.
%   b          - the observed vector of the size m x 1.
%   Aoper      - the linear operator.
%                It can be a matrix or a linear operator.
%                + If A is a matrix (real or complex), then its adjoint can 
%                  be computed automatically based on A. 
%                + If A is a linear operator, then the adjoint A' and the 
%                  number of variable nx must be provided.
%   param      - The optional parameters. It is a structure defined by user
%                See "lassoParamSettings" for more details.
%   varargin   - User defined inputs. It goes in pair (name, value)
%                   + Name = 'RegPar' - the regularization parameter.
%                     It is a positive scalar or nx-positive vector.
%                   + Name = 'NoiseLevel' - the noise level in the L1/L2con
%                     problem. It is a positive scalar.
%                   + Name = 'nx' - the number of variables.
%                   + Name = 'x0' - the initial point (nx - vector).
%                   + Name = 'AT' - the adjoint operator (function handle
%                     or matrix of the size nx x m).
%                   + Name = 'lbx' - the lower bound vector (nx - vector).
%                   + Name = 'ubx' - the upper bound vector (nx - vector).
%
% OUTPUTS:
%   xsol       - the final solution
%   output     - a structure containing all the performane information
%              including: 
%               + iter      = number of iterations
%               + time      = the cpu time
%               + rel_pfeas = the relative feasibility gap
%               + rel_schg  = the relative change of the solution
%               + fx_val    = the final objective value (if exists)
%               + msg       = the final status of the algorithm
%               + cntA      = number of matrix-vector multiplications Ax
%               + cntAt     = number of matrix-vector multiplications A'*y      
%
% SYNTAX:
%   For BP problem:
%     [xsol, output] = decoptLassoSolver('BP', A, b, param); 
%   or
%     [xsol, output] = decoptLassoSolver('BP', Aopr, b, param, 'AT', ATopr);
%
%   For L1/L2 problem:
%     [xsol, output] = decoptLassoSolver('L1/L2', A, b, param, 'RegPar', rho);
%
%   For L2/L1con problem:
%     [xsol, output] = decoptLassoSolver('L2/L1con', A, b, param, ...
%                                        'NoiseLevel', delta, 'lbx', 0);
%   and so on ...
% 
% INFORMATION:
%   By Quoc Tran Dinh, Laboratory for Informations and Inference Systems
%      (LIONS), EPFL, Lausanne, Switzerland.
%   Joint work with Volkan Cevher, LIONS, EPFL.
%   Date: 14.03.2014.
%   Last modified: 14.03.2014.
%   Contact: quoc.trandinh@epfl.ch
%   More information: http://homes.esat.kuleuven.be/~qtrandin/software.html
%
function [xsol, output] = decoptSolver(pstrType, A, b, param, varargin)

      
    % Check the inputs.
    if nargin < 4, param = []; end
    if nargin < 3, error('At least three inputs must be provided!'); end

    % Verify the optional parameters.
    param = decoptParamSettings(param); 
    param.PrintLength = 75;
    
    % Print the header file.
    decoptPrintMessages(param, 'header', param.PrintLength);

    % Preprocess the inputs.
    if param.Verbosity >= 1
        fprintf('%s\n', repmat('-', param.PrintLength, 1)); 
    end
    preprocess_time = tic;
    [pdata, opts] = decoptPreSolver(pstrType, A, b, varargin{:});

    % Print a message
    if param.Verbosity >= 0
        fprintf([' + DECOPT is solving ', opts.pstrType, ' problem ...\n']); 
    end

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
        kickPars.incrGam = 1.05;     % default: 1.05.
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