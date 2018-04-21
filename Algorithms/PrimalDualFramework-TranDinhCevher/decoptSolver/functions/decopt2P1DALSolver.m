% FUNCTION: [xsol, output] = decopt2P1DALSolver(Aopr, ATopr, b, ...
%                            proxOpers, lbx, ubx, x0, normAtA, param, ...
%                            opts, varargin)
%
% PURPOSE: An implementation of the decomposition algorithm for solving the
%          convex problem of the form:
%
%                       minimize_{x, r} f(x) + p(r)
%                       s.t.  A*x - r = b,
%
%          where f and p are two convex functions whose proximal operator
%          can be efficiently computed.
%
% ALGORITHM VARIANT:  Two Primal Steps and One Dual Step (2P1D) using
%                     Augmented Lagrangian smoothing technique.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 23.04.2014.
%    Contact: quoc.trandinh@epfl.ch
%
function [xsol, output] = decopt2P1DALSolver(Aopr, ATopr, b, proxOpers, ...
                          lbx, ubx, x0, normAtA, W, param, opts, varargin)

% Initialize the parameter gamma.
gamma      = 1e2*min(max(0.1/normAtA, opts.lbGamma), 1.0);  % default: 1e-6
if ~isempty(opts.gamma0), gamma = opts.gamma0; end
%opts.gamFact = 1;

% Initalize the parameters tau and beta.
tau        = 0.5*(sqrt(5) - 1);
beta       = tau^2/(gamma*(1 - tau));

% The Lipschitz constant of the dual function.
Lips       = normAtA;
LipsG      = 1;
mx         = length(b);
max_nrm_b1 = max(norm(b(:), 2), 1);

% The proximal operators.
fxProxOper = proxOpers{1};
prProxOper = proxOpers{2};
fxProxEval = proxOpers{3};
prProxEval = proxOpers{4};

% Generate the center points.
y_cen      = zeros(mx, 1);
r0         = -b;
cntA       = 0;
cntAt      = 0;
stopCrt    = 0;

% Compute the primal starting point xb.
[xb, rb, outb] = decoptCvxProbSolver(Aopr, ATopr, b, y_cen, ...
                 gamma, Lips, x0, r0, lbx, ubx, fxProxOper, ...
                 prProxOper, W, param, varargin{:});
xs         = xb;
rs         = rb;
subiter    = outb.iter;
cntA       = cntA  + outb.cntA;
cntAt      = cntAt + outb.cntAt;

% Compute the dual starting point yb.
Axb        = Aopr(xb);
Axh        = Axb;
cntA       = cntA + 1;
frstFeas   = Axb - rb - b;
yb         = (1.0/beta)*frstFeas;

% The main loop.
for iter = 1:param.MaxIters

    % STEP 1: Compute the primal variable xs and the residual rs.
    [xs, rs, outs] = decoptCvxProbSolver(Aopr, ATopr, b, yb, ...
                     gamma, Lips, xs, rs, lbx, ubx, fxProxOper, ...
                     prProxOper, W, param, varargin{:});
    subiter    = subiter + outs.iter;  
    cntA       = cntA    + outs.cntA;
    cntAt      = cntAt   + outs.cntAt;

    % STEP 2: Evaluate the objective value.
    if param.isFxEval
        fx_val = decoptFxEval(xb, rb, fxProxEval, prProxEval, xb, rb, ...
                              varargin{:});
    end

    % STEP 3: Update the primal variable xh and the residual rh.
    xh         = (1 - tau)*xb  + tau*xs;
    rh         = (1 - tau)*rb  + tau*rs;
    Axh_next   = Aopr(xh);
    cntA       = cntA + 1;
    
    % STEP 4: Update the smoothness parameter beta.
    beta       = (1 - tau)*beta;

    % STEP 5: Update the dual variable yh.
    frstFeas   = Axh_next - rh - b;
    yh         = (1/beta)*frstFeas;
    bh         = Axh_next;

    % STEP 6: Compute the primal variable xb_next and residual rb_next.
    [xb_next, rb_next, outb] = decoptCvxProbSolver(Aopr, ATopr, ...
                               bh, yh, (LipsG/beta), Lips, xb, rb, ...
                               lbx, ubx, fxProxOper, prProxOper, ...
                               W, param, varargin{:});
    subiter    = subiter + outb.iter;  
    cntA       = cntA    + outb.cntA;
    cntAt      = cntAt   + outb.cntAt;

    % STEP 7: Update the dual variable yb_next.
    yb_next    = (1 - tau)*yb + tau*yh;

    % STEP 8:  Computing the primal/dual feasibility and solution change.
    abs_pfeas  = norm(frstFeas(:), 2);
    abs_dfeas  = norm(Axh_next(:) - Axh(:), 2);
    abs_schg   = norm(xb_next(:) - xb(:), 2);
    rel_pfeas  = abs_pfeas/max_nrm_b1;
    rel_dfeas  = abs_dfeas/max(1.0, norm(Axh(:), 2));
    rel_schg   = abs_schg/max(1, norm(xb(:), 2));

    % STEP 9: Update the smoothness parameter gamma.
    decoptParamUpdate();
    %gamma_next = (1 - tau)^5*gamma_next;
    
    % STEP 10: Update the step size parameter tau.
    omeg       = gamma_next/gamma*tau^2;
    tau        = 0.5*(sqrt(omeg^2 + 4*omeg) - omeg);
    gamma      = gamma_next;

    % STEP 11: Print the iterations and save the history.
    decoptPrintIters();
    decoptSaveHistory();

    % STEP 12: Check the stopping criterion.
    if rel_pfeas <= param.RelTolFeas && rel_schg <= param.RelTolX
        stopCrt = 1;
        decoptTermination();
        break;
    end

    % STEP 13: Assign to the next iteration.
    xb    = xb_next;
    rb    = rb_next;
    yb    = yb_next;
    Axh   = Axh_next;
end

% Perform the final phase.
stopCrt   = 0;
decoptTermination();

% Get the final solution.
xsol             = xb_next;

% Get the final outputs.
output.stopCrt   = stopCrt;
output.iter      = iter;
output.rel_pfeas = rel_pfeas;
output.rel_dfeas = rel_dfeas;
output.rel_schg  = rel_schg;
output.fx_val    = fx_val;
output.cntA      = cntA;
output.cntAt     = cntAt;
output.auxi.rb   = rb_next;
output.auxi.yb   = yb_next;
output.auxi.xs   = xs;
output.auxi.rs   = rs;

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.