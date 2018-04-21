% FUNCTION: fx = decoptFxEval(x, r, fxProx, prProx, x0, r0, varargin)
% PURPOSE:  Evaluate the objective value of the problem:
%
%                        min_{x, r} f(x) + p(r)
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
function fx = decoptFxEval(x, r, fxProx, prProx, x0, r0, varargin)
    
    fx = fxProx(x, x0, varargin{:}) + prProx(r, r0, varargin{:});
    
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.