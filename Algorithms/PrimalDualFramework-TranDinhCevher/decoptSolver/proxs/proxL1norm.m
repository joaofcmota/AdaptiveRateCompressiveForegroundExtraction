% FUNCTION: zout = proxL1norm(zc, gamma)
% PURPOSE:  Compute the proximal operator of the l1-norm l(x) = ||x||_1.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 23.04.2014.
%    Contact: quoc.trandinh@epfl.ch
%
function zout = proxL1norm(zc, gamma)

    % Check the inputs.
    if nargin < 1,      gamma = 1; end
    if isempty(gamma),  gamma = 1; end
    
    % Perform the soft-thresholding operator.
    zout  = sign(zc).*max(abs(zc) - gamma, 0.0); 
   
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.