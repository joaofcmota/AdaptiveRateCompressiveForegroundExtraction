% FUNCTION: [pdata, opts] = decoptPreSolver(pstrType, A, b, varargin)
%
% PURPOSE:  Pre-process the inputs of the problem.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 14.03.2014.
%    Contact: quoc.trandinh@epfl.ch
%
function [pdata, opts] = decoptPreSolver(pstrType, A, b, varargin)

% Process problem type.
opts.pstrType = pstrType;
if strcmpi(pstrType, 'L1/L2'),     ts.P = 1;   ts.L1 = 1;  ts.Res = 1;  end
if strcmpi(pstrType, 'L1/L1'),     ts.P = 2;   ts.L1 = 1;  ts.Res = 2;  end
if strcmpi(pstrType, 'L1/sqrtL2'), ts.P = 3;   ts.L1 = 1;  ts.Res = 3;  end
if strcmpi(pstrType, 'L1/L2con'),  ts.P = 4;   ts.L1 = 1;  ts.Res = 4;  end
if strcmpi(pstrType, 'L2/L1con'),  ts.P = 5;   ts.L1 = 2;  ts.Res = 1;  end
if strcmpi(pstrType, 'BP'),        ts.P = 6;   ts.L1 = 1;  ts.Res = 4;  end
if strcmpi(pstrType, 'GBP'),       ts.P = 7;   ts.L1 = 3;  ts.Res = 4;  end
if strcmpi(pstrType, 'HingeL1'),   ts.P = 8;   ts.L1 = 1;  ts.Res = 5;  end
if strcmpi(pstrType, 'HingeL2'),   ts.P = 8;   ts.L1 = 4;  ts.Res = 5;  end
if strcmpi(pstrType, 'L1-L1'),     ts.P = 9;   ts.L1 = 5;  ts.Res = 4;  end
if strcmpi(pstrType, 'L1-L1NS'),   ts.P = 9;   ts.L1 = 5;  ts.Res = 4;  end
if strcmpi(pstrType, 'UserDef'),   ts.P = 100; ts.L1 = []; ts.Res = []; end

% Initialize the parameter pdata and opts.
opts.prox       = ts;
pdata.b         = b;
pdata.Aopr      = A;  pdata.ATopr  = []; 
pdata.nx        = []; pdata.x0     = [];
pdata.lbx       = []; pdata.ubx    = [];
pdata.label     = [];
pdata.varargin  = {};
opts.weights    = []; 
opts.delta      = 0;
opts.lambda     = 1;
groups          = [];
if opts.prox.P >= 4, opts.lambda  = 0; end

% Process the user-define parameters.
for idx = 1:2:(length(varargin)-1)
    switch upper(varargin{idx})
        case 'REGPAR'
            if isnumeric(varargin{idx+1}) 
                opts.weights = varargin{idx+1};
                if min(opts.weights(:)) < 0
                    error('The weighted vector must be positive!'); 
                end
            end
        case 'NX'
            if isnumeric(varargin{idx+1}) && ~isnan(varargin{idx+1}) ...
               && ~isinf(varargin{idx+1}) && varargin{idx+1} > 0
                pdata.nx = varargin{idx+1};
            end    
        case 'X0'
            if isnumeric(varargin{idx+1})
                pdata.x0 = varargin{idx+1};
            end
        case 'LBX'
            if isnumeric(varargin{idx+1})
                pdata.lbx = varargin{idx+1};
                if min(pdata.lbx(:)) > 0 || all(isinf(pdata.lbx))
                    error('The lower bound vector must be nonpositive!'); 
                end
            end
        case 'UBX'
            if isnumeric(varargin{idx+1})
                pdata.ubx = varargin{idx+1};
                if min(pdata.ubx(:)) < 0 || all(isinf(pdata.ubx))
                    error('The upper bound vector must be positive!'); 
                end
            end
        case 'NOISELEVEL'
            if isnumeric(varargin{idx+1}) && varargin{idx+1} > 0
                opts.delta = varargin{idx+1};
            end
        case 'AT'
            if isnumeric(varargin{idx+1}) || ...
               isa(varargin{idx+1}, 'function_handle')
               pdata.ATopr = varargin{idx+1};
            end    
        case 'GROUPS'
            if isstruct(varargin{idx+1})
               groups = varargin{idx+1};
            else
               error('The groups should be a structure!'); 
            end
        case 'LABEL'
            if isnumeric(varargin{idx+1})
                pdata.label = varargin{idx+1};
            end
        case 'PROX'    
            if iscell(varargin{idx+1}) && length(varargin{idx+1}) == 4
                pdata.proxOpers = varargin{idx+1};
            else
                error('Prox-functions are not correctly provided!');
            end
        case 'GAMMAFACTOR'
            if isnumeric(varargin{idx+1}) && varargin{idx+1} >= 1
                opts.incrGam = varargin{idx+1};
            end
            
        case 'GAMMA0'
            if isnumeric(varargin{idx+1}) && varargin{idx+1} > 0
                opts.gamma0 = varargin{idx+1};
            end
     end
end

% Check the linear operator A and its adjoint.
if isnumeric(A) 
    pdata.Aopr  = @(x)(A*x); 
    pdata.ATopr = @(y)(A'*y); 
    pdata.nx    = size(A, 2);
end

% If Ax is a function handle 
if isa(A, 'function_handle')
    pdata.Aopr = A;
    if isempty(pdata.ATopr)
        error('The adjoint operator must be provided!');
    end
    if isempty(pdata.nx)
        error('The number of variables must be provided!');
    end
end

% Now, check the size consistancy.
if isempty(pdata.nx) || pdata.nx < 1
    error('The number of variable is incorrect'); 
end
if isempty(pdata.x0) && ~isempty(pdata.nx)
    pdata.x0 = zeros(pdata.nx,1); 
end
if ~isempty(pdata.x0) && length(pdata.x0) ~= pdata.nx
    error('The initial point is not correct!');
end
if ~isempty(pdata.ubx) && (length(pdata.ubx) ~= pdata.nx || ...
    all(isinf(pdata.ubx)))
    error('The upper bound is not correct!');
end
if ~isempty(opts.weights) && length(opts.weights) > 1 ...
    && ( length(opts.weights) ~= pdata.nx && opts.prox.P ~= 7 )
    error('The weighted vector is not correct!');
end
if ~isempty(opts.delta) && ~isnumeric(opts.delta) && opts.delta <= 0
    error('The selector level delta is not correct');
end
if isempty(opts.weights), opts.weights = 1.0; end

% Make sure that every inputs are vectors.
pdata.x0     = pdata.x0(:);
pdata.lbx    = pdata.lbx(:);
pdata.ubx    = pdata.ubx(:);
opts.weights = opts.weights(:); 

% Check if A is orthonormal.
xtmp = randn(pdata.nx, 1);
ytmp = pdata.ATopr(pdata.Aopr(xtmp));
err  = norm(ytmp - xtmp, 2)/max(1, norm(xtmp, 2));
if err <= 1e-9, opts.isOrtA = true; else opts.isOrtA = false; end

% Check if prox-operators are provided.
if opts.prox.P == 100 && isempty(pdata.proxOpers)
    error('The prox functions are not provided!');
end

% Preprocess the groups
if ~isempty(groups)
    pdata.ng      = groups.ng;
    pdata.groups  = groups.indices;
    pdata.mgroups = sparse(groups.ng, pdata.nx);
    for i=1:groups.ng
        pdata.mgroups(i, groups.indices == i) = 1;
    end
else
    pdata.ng = 0; pdata.groups = []; pdata.mgroups = [];
end
    
% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.
