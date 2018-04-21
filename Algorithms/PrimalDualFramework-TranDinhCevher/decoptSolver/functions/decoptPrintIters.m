% FUNCTION: decoptPrintIters()
% PURPOSE:  Print the iteration if required.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
%% FUNCTION: decoptPrintIters()

% Not output is printed.
if param.Verbosity <= 1, return; end

% Print the header.
if mod(iter, 10*param.PrintStep) == 1 || iter == 1
    fprintf('%s\n', repmat('-', 1, param.PrintLength)); 
    fprintf([' Iter| RelPGap | RelDGap | RelSchg |', ...
             '  Beta  | Gamma  |   Tau  |    F(x)\n']);
    fprintf('%s\n', repmat('-', 1, param.PrintLength));
end

% Print the values: iter -> rel_pfeas -> rel_dfeas -> rel_schg 
%                   beta -> gamma -> tau -> fx_val
if mod(iter, param.PrintStep) == 0 || iter == 1
    fprintf(['%5d| %3.2e| %3.2e| %3.2e| %3.1e| %3.1e|', ...
             ' %3.1e| %3.5e \n'], ...
             iter, rel_pfeas, rel_dfeas, rel_schg, ...
             beta, gamma, tau, fx_val);
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.