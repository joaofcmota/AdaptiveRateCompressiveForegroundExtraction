% FUNCTION: [xsol, output] = decopt2P1DsSolver(Aopr, ATopr, b, ...
%                            proxOpers, lbx, ubx, x0, normAtA, W, ...
%                            sigma_f, param, opts, varargin)
%
% PURPOSE: An implementation of the decomposition algorithm for solving the
%          convex problem of the form:
%
%                       minimize_{x, r} f(x) + p(r)
%                       s.t.  A*x - r = b,
%
%          where f and p are two convex functions whose proximal operator
%          can be efficiently computed. Here, we also assume that f is 
%          strongly convex with a parameter sigma_f > 0.
%
% ALGORITHM VARIANT:  Two Primal Steps and One Dual Step (2P1D) using
%                     Bregman distance smoothing technique.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 23.04.2014.
%    Contact: quoc.trandinh@epfl.ch
%   
function [xsol, output] = decopt2P1DsSolver(Aopr, ATopr, b, proxOpers, ...
                          lbx, ubx, x0, normAtA, W, sigma_f, ...
                          param, opts, varargin)

% This is the Lipschitz constant of the dual function g(y).
lipsG    = normAtA;
dualLips = lipsG + 1;

% Initialize the parameters tau, gamma and beta.
tau      = 0.5*(sqrt(5) - 1);
gamma    = sigma_f;
beta     = dualLips*tau^2/(gamma*(1 - tau));

% Define the proximal operator of f and p.
fxProxOpr  = proxOpers{1};
prProxOpr  = proxOpers{2};
fxProxEval = proxOpers{3};
prProxEval = proxOpers{4};

% Initialize the outputs.
output   = []; 
cntA     = 0;
cntAt    = 0; 
fx_val   = nan;
norm_b   = norm(b, 2);
stopCrt  = 0;

% Define the center points.
x_cen    = x0;
y_cen    = zeros(size(b));
Ax_cen   = Aopr(x_cen);
cntA     = cntA + 1;

% Compute the primal residual rb.
igamma   = 1/gamma;
rg       = Ax_cen - b + igamma*y_cen;
rs       = prProxOpr(rg, igamma, rg, varargin{:});
rb       = rs;

% Define the primal point xs and rs.
if param.adaptStepSize
    ATyb  = gamma*ATopr((Ax_cen - rs - b) + igamma*y_cen);
    fg    = ATyb;
    cntAt = cntAt + 1;
    Afg   = Aopr(fg);
    cntA  = cntA  + 1;
    alpha = (fg'*fg)/(gamma*(Afg'*Afg) + eps);
    xg    = x_cen - alpha*fg;
else
    iLips = 1.0/lipsG;
    igamL = 1.0/(gamma*lipsG);
    ATyb  = ATopr(igamL*y_cen + iLips*(Ax_cen - rs - b));
    cntAt = cntAt + 1;
    xg    = x_cen - ATyb;
    alpha = igamL;
end
xs       = fxProxOpr(xg, alpha.*W, x_cen, varargin{:});
if ~isempty(lbx), xs = max(xs, lbx); end
if ~isempty(ubx), xs = min(xs, ubx); end
xb       = xs;

% Applying the linear operator.
Axb      = Aopr(xb);
cntA     = cntA + 1;

% Update the dual variable yb.
frstFeas = Axb - rb - b;
yb       = y_cen + (1.0/beta)*frstFeas;

% The main loop of the algorithm.
for iter = 1:param.MaxIters
    
    % STEP 1: Evaluate the objective value.
    if param.isFxEval
        fx_val = decoptFxEval(xb, rb, fxProxEval, prProxEval, xb, rb, ...
                 varargin{:});
    end
    
    % STEP 2: Compute the primal point xs and rs - First step.
    igamma     = 1.0/gamma;
    rg        = rs + igamma*yb;
    rs_next   = prProxOpr(rg, igamma, rb, varargin{:});
    
    % STEP 3: Compute the primal solution xs.
    xg        = xb - igamma*ATyb;
    alpha     = igamma;
    xs_next   = fxProxOpr(xg, alpha.*W, xb, varargin{:});
    if ~isempty(lbx), xs_next = max(xs_next, lbx); end
    if ~isempty(ubx), xs_next = min(xs_next, ubx); end

    % STEP 4: Applying the linear operator A.
    Axs_next  = Aopr(xs_next);
    cntA      = cntA + 1;
    
    % STEP 5: Update xh.
    rh        = (1 - tau)*rb  + tau*rs_next;
    xh        = (1 - tau)*xb  + tau*xs_next;
    Axh       = (1 - tau)*Axb + tau*Axs_next;
    
    % STEP 6: Update the parameter beta.
    beta      = (1 - tau)*beta;
    
    % STEP 7: Update the dual variable yh.
    frstFeas  = Axh - rh - b;
    yh        = y_cen + (1.0/beta)*frstFeas;
    
    
    % STEP 8: Compute the primal point xs and rs - Second step.
    igamma2   = beta/lipsG;
    ATyh      = ATopr(yh);
    cntAt     = cntAt + 1;
    xgh       = xh - igamma2*ATyh;
    alp2      = igamma2;
    xb_next   = fxProxOpr(xgh, alp2.*W, xh, varargin{:});
    if ~isempty(lbx), xb_next = max(xb_next, lbx); end
    if ~isempty(ubx), xb_next = min(xb_next, ubx); end
    
    % STEP 9: Applying the linear operator A.        
    Axb_next  = Aopr(xb_next);
    cntA      = cntA + 1;
    
    % STEP 10: Compute the residual term - Second step.
    %rgh       = Axb_next - b + igamma2*yh;
    rgh       = rh + igamma2*yh;
    rb_next   = prProxOpr(rgh, igamma2, rh, varargin{:});
    sndFeas   = Axb_next - rb_next - b;
    
    % STEP 11: Update yb_next.
    yb_next   = (1 - tau)*yb + tau*yh;
    ATyb_next = (1 - tau)*ATyb + tau*ATyh;
    
    % STEP 12: Evaluate the absolute primal and dual feasibility.
    abs_schg  = norm(xb_next(:) - xb(:), 2);
    abs_pfeas = norm(sndFeas(:), 2);
    abs_dfeas = gamma*norm(Axb_next(:) - Axb(:), 2);
    
    % STEP 13: Evaluate the relative primal and dual feasibility.
    rel_pfeas = abs_pfeas/max(1.0, norm_b);
    rel_dfeas = abs_dfeas/max(1.0, norm(rb(:), 2));
    rel_schg  = abs_schg/max(norm(xb, 2), 1.0);
    
    % STEP 14: Update the smoothness parameter gamma.
    decoptParamUpdate();
    gamma_next = (1 - tau)^(1/3)*gamma_next;
    
    % STEP 15: Print the iteration and save the history if required.
    decoptPrintIters();
    decoptSaveHistory();
    
    % STEP 16: Check the stopping criterion.
    if rel_pfeas <= param.RelTolFeas && rel_schg <= param.RelTolX
        stopCrt  = 1;
        decoptTermination();
        break;
    end
    
    % STEP 17: Update the parameter tau.
    omeg  = gamma_next/gamma*tau^2;
    tau   = 0.5*(sqrt(omeg^2 + 4*omeg) - omeg);
    gamma   = gamma_next;
    
    % STEP 18: Assign to the next iteration.
    xs    = xs_next;
    rs    = rs_next;
    xb    = xb_next;
    yb    = yb_next;
    rb    = rb_next;
    Axb   = Axb_next;
    ATyb  = ATyb_next;
end
% End of the main loop.

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