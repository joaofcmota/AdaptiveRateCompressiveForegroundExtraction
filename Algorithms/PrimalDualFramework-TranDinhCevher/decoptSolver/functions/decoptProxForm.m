% FUNCTION: proxOpers = decoptProxForm(data, opts)
%
% PURPOSE:  Process the prox operators.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 14.03.2014.
%    Contact: quoc.trandinh@epfl.ch
%
function proxOpers = decoptProxForm(data, opts)

% The weighted vector and the radius.
W       = opts.weights;
radius  = opts.delta;
groups  = data.groups;
mgroups = data.mgroups;
label   = data.label;
w       = data.w;
beta    = data.beta;

% -------------------------------------------
% The proximal operator for the x-term: f(x).
% -------------------------------------------
% The l1-norm
if opts.prox.L1 == 1
    proxOpers{1} = @(x, gamma, varargin)(proxL1norm(x, gamma));
    proxOpers{3} = @(x, varargin)( norm(W.*x(:), 1) );
end


% The projection on L1-norm ball.
if opts.prox.L1 == 2
    proxOpers{1} = @(x, gamma, varargin)(projL1norm(x, radius));
    proxOpers{3} = @(x, varargin)(0);
end

% The proximal operator of group lasso.
if opts.prox.L1 == 3
    proxOpers{1} = @(x, gamma, varargin)(proxGroupL1L2norm(x, gamma, ...
                     groups, mgroups, 1.0));
    proxOpers{3} = @(x, varargin)(sum(W.*sqrt(mgroups*x(:).^2)));
end

% The proximal operator of the square l2-norm.
if opts.prox.L1 == 4
    proxOpers{1} = @(x, gamma, varargin)(proxSquareL2norm(x, gamma));
    proxOpers{3} = @(x, varargin)(0.5*sum(W.*x.^2));
end

% The L1L1-norm
if opts.prox.L1 == 5
    proxOpers{1} = @(x, gamma, varargin)(proxL1L1norm(x, w, beta, gamma));             % CREATE THIS FUNCTION
    proxOpers{3} = @(x, varargin)( norm(x(:), 1) + beta*norm(x(:) - w(:), 1));
end

% ----------------------------
% For the residual term: p(r).
% ----------------------------
% The square L2-norm.
if opts.prox.Res == 1 
    proxOpers{2} = @(x, gamma, varargin)(proxSquareL2norm(x, gamma));
    proxOpers{4} = @(x, varargin)(0.5*norm(x(:), 2)^2);
end

% The L1-norm.
if opts.prox.Res == 2 
    proxOpers{2} = @(x, gamma, varargin)(proxL1norm(x, gamma));
    proxOpers{4} = @(x, varargin)(norm(x(:), 1));
end

% The sqrt L2-norm.
if opts.prox.Res == 3
    proxOpers{2} = @(x, gamma, varargin)(proxL2norm(x, gamma));
    proxOpers{4} = @(x, varargin)(norm(x(:), 2));
end

% The projection on L2-norm ball.
if opts.prox.Res == 4
    proxOpers{2} = @(x, gamma, varargin)(projL2norm(x, radius));
    proxOpers{4} = @(x, varargin)(0);
end

% The Hinge-loss proximal operators.
if opts.prox.Res == 5
    proxOpers{2} = @(x, gamma, varargin)(proxHingeLoss(x, gamma, label, label.^2));
    proxOpers{4} = @(x, varargin)(sum(max(1 - label.*x, 0)));
end

% User-define prox functions.
if opts.prox.P == 100 && iscell(data.proxOpers) && length(data.proxOpers) == 4
    for ii=1:4
        if isa(data.proxOpers{ii}, 'function_handle')
            proxOpers{ii} = data.proxOpers{ii};
        else
            error('The prox function is not valid!');
        end
    end
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.
