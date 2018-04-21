% FUNCTION: [xsol, rsol, output] = decoptCvxProbSolver(Aoper, AToper, b, ...
%                                  ys, gamma, Lips, x0, r0, lbx, ubx, ...
%                                  fxProxOper, prProxOper, W, ...
%                                  param, varargin)
%
% PURPOSE: Solve the convex subproblem of the form:
%
%  min_{x, r} f(x) + p(r) + y'*(A*x - r - b) + 0.5*gamma*|A*x - r - b|^2_2.
%
%  where gamma > 0 is the penalty parameter and y is the multipliers.
% 
% METHOD: An implementation of the FISTA algorithm for solving this problem
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 23.04.2014.
%    Contact: quoc.trandinh@epfl.ch 
%
function [xsol, rsol, output] = decoptCvxProbSolver(Aoper, AToper, b, ...
                                ys, gamma, Lips, x0, r0, lbx, ubx, ...
                                fxProxOper, prProxOper, W, param, varargin)
                                       
% Initialize the control parametes.
cntA  = 0;
cntAt = 0;
time1 = tic;

% Intermediate parameters.
igamma   = 1.0./gamma;
igam_ys  = igamma.*ys;
iLips    = 1/Lips;
igamLips = W./(gamma*Lips);

% Initialize the variables.
x1_cur = x0;
x2_cur = x1_cur;
r1_cur = r0;
r2_cur = r1_cur;
t_cur  = 1;

% The main loop.
for iter = 1:param.InnerMaxIters

    % Applying the linear operator.
    Ax2_cur  = Aoper(x2_cur);
    cntA     = cntA + 1;

    % Applying the proximal operator on the r-term.
    r_cen    = Ax2_cur - b + igam_ys;
    r1_next  = prProxOper(r_cen, igamma, r1_cur, varargin{:});

    % Applying the proximal operator on the x-term.
    x_ss     = (Ax2_cur - r1_next - b) + igam_ys;
    %x_ss     = (Ax2_cur - r2_cur - b) + igam_ys;
    ATy      = AToper(x_ss);
    cntAt    = cntAt + 1;
    x_cen    = x2_cur - iLips*ATy;
    x1_next  = fxProxOper(x_cen, igamLips, x1_cur, varargin{:});

    % Project onto the box constraints.
    if ~isempty(lbx), x1_next = max(x1_next, lbx); end
    if ~isempty(ubx), x1_next = min(x1_next, ubx); end

    % Update the parameter t_cur.
    t_next   = 0.5*(sqrt(t_cur^2 + 4) + 1);

    % Update the variable z_cur.
    step_cur = (t_cur - 1)/t_next;
    x2_next  = x1_next + step_cur*(x1_next - x1_cur);
    r2_next  = r1_next + step_cur*(r1_next - r1_cur);

    % Compute the solution changes.
    abs_schg = norm([x1_next(:) - x1_cur(:); r1_next(:) - r1_cur(:)], 2);
    rel_schg = abs_schg/max(1, norm([x1_cur(:); r1_cur(:)], 2));

    % Check stopping criterion.
    if rel_schg <= param.InnerRelTol 
        output.msg = 'Convergence achieved'; 
        break;
    end

    % Perform fixed restart.
    if mod(iter, 50) == 0
        x2_next = x1_next;
        r2_next = r1_next;
        t_next  = 1;
    end

    % Move to the next iteration.
    x1_cur = x1_next;
    x2_cur = x2_next;
    r1_cur = r1_next;
    r2_cur = r2_next;
    t_cur  = t_next;
end

% Perform the final phase.
if iter >= param.InnerMaxIters
    output.msg = 'Exceed the maximum number of iterations';
end

% Get the final solution.
xsol            = x1_next;
rsol            = r1_next;
output.iter     = iter;
output.abs_schg = abs_schg;
output.rel_schg = rel_schg;
output.cntA     = cntA;
output.cntAt    = cntAt;
output.time     = toc(time1);

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.