% FUNCTION: zout = proxHingeLoss(zc, gamma, yf, yf2)
% PURPOSE:  Compute the proximal operator of the Hinge loss function:
%
%                  H(t, yf) := max{1 - yf*t, 0) = [1 - yf*t]_{+}
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 23.04.2014.
%    Contact: quoc.trandinh@epfl.ch
%
function zout = proxHingeLoss(zc, gamma, yf, yf2)

    % Check the inputs.
    if nargin < 3,   error('At least three inputs must be provided!'); end
    if nargin < 4,   yf2 = yf.^2; end
    if isempty(yf2), yf2 = yf.^2; end

    % Compute the threshold value.
    beta = 0.5*gamma.*yf2;
    
    % Shift the center points.
    zcb  = (1.0 - yf.*zc) - beta;
    
    % Perform the soft-thresholding operator.
    sss  = sign(zcb).*max(abs(zcb) - beta, 0.0); 
    
    % Recover the solution zout.
    zout = (1.0 - sss)./yf;
    
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.