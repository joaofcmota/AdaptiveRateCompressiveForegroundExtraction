% FUNCTION: xout = proxGroupL1L2norm(ss, gamma, groups, mgroups, W)
% PURPOSE:  Perform the proximal operator of the L1-L2-group norm:
%                    min_x g(x) + 1/(2*gamma)*|x - s|_2^2,
%           where g(x) = sum_i=1^ng gamma(i)*|x(i)|_2.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
function xout = proxGroupL1L2norm(ss, gamma, groups, mgroups, W)
    
    if nargin < 4, error('At least four inputs must be provided!'); end
    if isempty(gamma), gamma = 1; end
    if nargin > 4, gamma = gamma.*W; end

    ss_group = max(0, 1 - gamma./sqrt(mgroups*ss.^2));
    xout     = ss_group(groups).*ss;
    
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.