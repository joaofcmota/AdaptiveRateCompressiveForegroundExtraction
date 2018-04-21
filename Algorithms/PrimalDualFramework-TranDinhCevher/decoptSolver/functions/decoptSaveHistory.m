% FUNCTION: decoptSaveHistory()
% PURPOSE:  Save the history w.r.t. the iterations if requires.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 23.04.2014
%    Contact: quoc.trandinh@epfl.ch
%
%% FUNCTION: decoptSaveHistory()

% Level 1: Basic information: relative error, residual and objective values.
if param.saveHistMode > 0
    output.hist.rel_schg(  iter, 1) = rel_schg;
    output.hist.rel_pfeas( iter, 1) = rel_pfeas;
    output.hist.rel_dfeas( iter, 1) = rel_dfeas;
    output.hist.fx_val(    iter, 1) = fx_val;
end

% Level 2: Absolute error and residual.
if param.saveHistMode > 1
    output.hist.abs_schg(  iter, 1) = abs_schg;
    output.hist.abs_pfeas( iter, 1) = abs_pfeas;
    output.hist.abs_dfeas( iter, 1) = abs_dfeas;
end

% Level 3: Parameters.
if param.saveHistMode > 2
   output.hist.beta(iter, 1)  = beta;
   output.hist.gamma(iter, 1) = gamma;
   output.hist.tau( iter, 1)  = tau;
end

% Level 3: Save the information for the inner loop.
if param.saveHistMode > 3
    if exist('outs')
        if isfield(outs, 'iter')
            output.hist.subit(iter, 1)   = outs.iter;
        end
        if isfield(outs, 'time')
            output.hist.subtime(iter, 1) = outs.time;
        end
        if isfield(outs, 'rel_schg')
            output.hist.suberr(iter, 1)  = outs.rel_schg;
        end
    end
end

% Level 5: Iterative vectors.
if param.saveHistMode > 4
   output.hist.xb{iter}       = xb;
   output.hist.rb{iter}       = rb;
   output.hist.yb{iter}       = yb;
   output.hist.xs{iter}       = xs;
   output.hist.rs{iter}       = rs;
end
       
% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.