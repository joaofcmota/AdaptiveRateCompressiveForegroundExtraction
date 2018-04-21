% FUNCTION: decoptPrintMessages(param, mode, varargin) 
% PURPOSE:  Print the output information if require.
% 
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
function decoptPrintMessages(param, mode, varargin)

% Not output is printed.
if param.Verbosity <= 0, return; end

% Print the header.
if strcmpi(mode, 'header') 
    if nargin > 2, n_sp = varargin{1}; else n_sp = 72; end
    fprintf('%s\n\n', repmat('*', 1, n_sp));
    fprintf('')            
    fprintf('       -----------------------------------------------------\n');
    fprintf('       A PRIMAL-DUAL DECOMPOSITION FRAMEWORK FOR CONSTRAINED\n');
    fprintf('                       CONVEX OPTIMIZATION \n');
    fprintf('       -----------------------------------------------------\n\n');
    fprintf(' This is a MATLAB software package for solving l1-norm convex\n'); 
    fprintf(' minimization problems of a linear meassurment operator.\n');
    fprintf(' It consists of several algorithmic variants for solving \n');
    fprintf(' the following different problems: L1/L2, L1/L1, L1/sqrtL2\n');
    fprintf(' L1/L2con, L2/L1con, Basis Pursuit, Group-Sparse LASSO, Binary SVM\n');
    fprintf(' with Hingle Loss ...\n');
    fprintf(' For more information, we refer to our website: http://lion.epfl.ch.\n\n');
    
    fprintf(' CONTACT INFORMATION:\n');
    fprintf('  By Quoc Tran-Dinh\n');
    fprintf('    Laboratory for Information and Inference Systems (LIONS)\n');
    fprintf('    EPFL, Lausanne, Switzerland.\n');
    fprintf('  Joint work with Volkan Cevher (LIONS).\n');
    fprintf('  Contact: quoc.trandinh@epfl.ch\n');
    fprintf('  Date: 14.03.2014.\n');
    fprintf('  Last modified: 24.04.2014.\n\n');
end

% Print the final outputs.
if strcmpi(mode, 'output')
    if nargin > 2, fout = varargin{1}; end
    if nargin > 3, n_sp = varargin{2}; else n_sp = 60; end
    n_sp2 = round(max(n_sp/2 - 8, 0)); 
    fprintf('%s', repmat('*', 1, n_sp2));
    fprintf('THE FINAL OUTPUT');
    fprintf('%s\n', repmat('*', 1, n_sp2));

    if isfield(fout, 'alg')
        fprintf(' + Algorithm: %s\n', fout.alg);
    end
    if isfield(fout, 'status')
        fprintf(' + Status:    %s\n', fout.status);
    end
    if isfield(fout, 'iter')
        fprintf(' + Number of total iterations:       %d\n', fout.iter);
    end
    if isfield(fout, 'cntA')
        fprintf(' + The number of linear operations:  %d\n', fout.cntA);
    end
    if isfield(fout, 'cntAt')
        fprintf(' + The number of adjoint operations: %d\n', fout.cntAt);
    end
    if isfield(fout, 'total_time')
        fprintf(' + CPU time (s):                     %0.5f\n', fout.total_time);
        fprintf('   - Including (s): [%0.3f(preprocess), %0.3f(solve)]\n', ...
                    fout.pre_time, fout.solve_time);
    end
    if isfield(fout, 'fx_val') && ~isnan(fout.fx_val) && ~isinf(fout.fx_val)
        fprintf(' + The final objective value:           %0.12f\n', fout.fx_val);
    end
    if isfield(fout, 'rel_pfeas')
        fprintf(' + The relative primal feasibility gap: %0.12f\n', fout.rel_pfeas);
    end
    if isfield(fout, 'rel_dfeas')
        fprintf(' + The relative dual feasibility gap:   %0.12f\n', fout.rel_dfeas);
    end
    if isfield(fout, 'rel_schg')
        fprintf(' + The relative solution change:        %0.12f\n', fout.rel_schg);
    end
    fprintf('%s\n', repmat('*', 1, n_sp));
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.