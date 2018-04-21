% FUNCTION: decoptParamUpdate()
% PURPOSE:  Update the first smoothness parameter gamma.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
%% FUNCTION: decoptParamUpdate()

if abs_pfeas >= opts.gamFact*abs_dfeas
    gamma_next = gamma*opts.incrGam;
elseif abs_dfeas >= opts.gamFact*abs_pfeas
    gamma_next = gamma/opts.decrGam;
else
    gamma_next = gamma;
end

gamma_next = min(max(gamma_next, opts.minGam), opts.maxGam);

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.