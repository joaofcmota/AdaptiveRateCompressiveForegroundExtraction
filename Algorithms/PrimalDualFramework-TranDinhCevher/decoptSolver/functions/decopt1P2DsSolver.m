% FUNCTION: [fsol, output] = decopt1P2DsSolver(Aopr, ATopr, b, ...
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
%          can be efficiently computed. The function f is assumed to be 
%          strongly convex with parameter sigma_f.
%
% ALGORITHM VARIANT:  One Primal Step and Two Dual Steps (1P2D) using
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
function [xsol, output] = decopt1P2DsSolver(Aopr, ATopr, b, ...
                          proxOpers, lbx, ubx, x0, normAtA, W, ...
                          sigma_f, param, opts, varargin)
                      
% This is the Lipschitz constant of the dual function g(y).
lipsG    = normAtA;
dualLips = lipsG + 1;

% Initialize the parameters tau, gamma and beta.
tau      = 0.5*(sqrt(5) - 1);
gamma    = sigma_f;
%if ~isempty(opts.gamma0), gamma = opts.gamma0; end
beta     = dualLips*tau^2/(gamma*(1 - tau));

% The parameter for updating dual step.
dFactor  = 0.5*(sqrt(5) + 1);

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
    fg    = gamma*ATopr((Ax_cen - rs - b) + (1/gamma)*y_cen);
    cntAt = cntAt + 1;
    Afg   = Aopr(fg);
    cntA  = cntA  + 1;
    alpha = (fg'*fg)/(gamma*(Afg'*Afg) + eps);
    xg    = x_cen - alpha*fg;
else
    iLips = 1.0/lipsG;
    igamL = 1.0/(gamma*lipsG);
    xg    = x_cen - ATopr(igamL*y_cen + iLips*(Ax_cen - rs - b));
    cntAt = cntAt + 1;
    alpha = igamL;
end
xs       = fxProxOpr(xg, alpha.*W, x_cen, varargin{:});
if ~isempty(lbx), xs = max(xs, lbx); end
if ~isempty(ubx), xs = min(xs, ubx); end
xb       = xs;

% Applying the linear operator.
Axb      = Aopr(xb);
cntA     = cntA + 1;
Axs      = Axb;

% Update the dual variable yb.
frstFeas = Axb - rb - b;
yb       = y_cen + (1.0/beta)*frstFeas;

% The main loop of the algorithm.
for iter = 1:param.MaxIters
    
    % STEP 1: Evaluate the first dual variable ys.
    frstFeas  = Axb - rb - b;
    ybs       = y_cen + (1.0/beta)*frstFeas;
    
    % STEP 2: Compute the initermediate point yh.
    yh        = (1.0 - tau)*yb + tau*ybs;
    
    % STEP 3: Evaluate the objective value.
    if param.isFxEval
        fx_val = decoptFxEval(xb, rb, fxProxEval, prProxEval, ...
                              xs, rs, varargin{:});
    end
    
    % STEP 4: Compute the primal residual rs.
    igamma    = 1.0/gamma;
    rg        = Axs - b + igamma*yh;
    rs_next   = prProxOpr(rg, igamma, rs, varargin{:});

    % STEP 5: Compute the primal solution xs.
    if param.adaptStepSize
        fg    = gamma*ATopr((Axs - rs_next - b) + (1/gamma)*yh);
        cntAt = cntAt + 1;
        Afg   = Aopr(fg);
        cntA  = cntA  + 1;
        alpha = (fg'*fg)/(gamma*(Afg'*Afg) + eps);
        xg    = xs - alpha*fg;
    else
        iLips = 1.0/lipsG;
        igamL = 1.0/(gamma*lipsG);
        xg    = xs - ATopr(igamL*yh + iLips*(Axs - rs_next - b));
        cntAt = cntAt + 1;
        alpha = igamL;
    end
    xs_next   = fxProxOpr(xg, alpha.*W, xs, varargin{:});
    if ~isempty(lbx), xs_next = max(xs_next, lbx); end
    if ~isempty(ubx), xs_next = min(xs_next, ubx); end

    Axs_next  = Aopr(xs_next);
    cntA      = cntA + 1;
    
    % STEP 6: Update the dual variable.
    step      = dFactor*gamma;
    sndFeas   = Axs_next - rs_next - b;
    yb_next   = yh + step*sndFeas;
    %yb_next   = yb + step*sndFeas;
    
   
    % STEP 7: Update xb_next and rb_next.
    xb_next   = (1 - tau)*xb + tau*xs_next;
    rb_next   = (1 - tau)*rb + tau*rs_next;
    
    % STEP 8: Evaluate the primal and dual feasibility.
    ATd       = ATopr(rs_next - rs);
    cntAt     = cntAt + 1;
    abs_schg  = norm(xb_next(:) - xb(:), 2);
    abs_pfeas = norm(sndFeas(:), 2);
    abs_dfeas = gamma*norm(ATd(:), 2);
    
    rel_pfeas = abs_pfeas/max(1.0, norm_b);
    rel_dfeas = abs_dfeas/max(1.0, norm(rs(:), 2));
    rel_schg  = abs_schg/max(norm(xb, 2), 1.0);
    
    % STEP 9a: Update the smoothness parameter beta.
    beta      = (1 - tau)*beta;
    
    % STEP 9b: Update the smoothness parameter gamma.
    %decoptParamUpdate();
    
    % STEP 10: Print out the iterations and save the history.
    decoptPrintIters();
    decoptSaveHistory();
    
    % STEP 11: Check the stopping criterion.
    if rel_pfeas <= param.RelTolFeas && rel_schg <= param.RelTolX
        stopCrt = 1;
        decoptTermination();
        break;
    end
    
    % STEP 12: Update the parameter tau.
    gamma_next  = gamma;
    omeg        = gamma_next/gamma*tau^2;
    tau         = 0.5*(sqrt(omeg^2 + 4*omeg) - omeg);
    gamma       = gamma_next;
    
    % STEP 13: Assign to the next iteration.
    rs    = rs_next;
    xs    = xs_next;
    yb    = yb_next;
    xb    = xb_next;
    rb    = rb_next;
    Axs   = Axs_next;
    Axb   = (1.0 - tau)*Axb + tau*Axs;
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