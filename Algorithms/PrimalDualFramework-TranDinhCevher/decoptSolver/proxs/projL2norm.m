% FUNCTION: xout = projL2norm(ss, gamma)
% PURPOSE:  Perform the projection operator on the l2 norm ball.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
function xout = projL2norm(ss, gamma)

    if nargin < 2,     gamma = 1; end
    if isempty(gamma), gamma = 1; end
    
    norm_ss = norm(ss(:), 2);
    xout    = gamma*min(1, 1/norm_ss)*ss;
    
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.