% FUNCTION: decoptTermination()
% PURPOSE:  Perform the termination procedure.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 23.04.2014
%    Contact: quoc.trandinh@epfl.ch
%
%% FUNCTION: decoptTermination()

% CASE 1: Convergence achieved if the feasible gap & solution changes are 
%         below the given tolerance.
if stopCrt == 1
    output.status = 'Convergence achieved';
    output.msg    = ['Feasibility gap and diffences of solutions', ...
                     ' are below the given tolerance and the ',...
                     'search direction does not change signficantly'];

% CASE 2: Convergence achieved if the feasible gap & difference of the 
%         objective values are below the given tolerance.
elseif stopCrt == 2
    output.status = 'Convergence achieved';
    output.msg    = ['Feasibility gap is below the given tolerance',...
                     ' and the objective does not change signficantly'];
              
% CASE 3: Reach the maximum number of iterations.           
elseif stopCrt == 0 && iter >= param.MaxIters
    output.status = 'Have not converged yet';
    output.msg    = ['Exceed the maximum number of iterations. ', ...
                     'Increase MaxIters if required!'];
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.