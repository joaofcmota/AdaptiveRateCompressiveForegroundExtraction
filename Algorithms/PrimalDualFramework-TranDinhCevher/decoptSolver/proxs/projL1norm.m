% FUNCTION: xout = projL1norm(x, delta)
% PURPOSE:  Perform the projection on the L1-ball.
%
% INFORMATION:
%    By Quoc Tran-Dinh, Laboratory for Informations and Inference Systems
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 14.03.2014
%    Last modified: 14.03.2014
%    Contact: quoc.trandinh@epfl.ch
%
function xout = projL1norm(x, delta)

    if nargin < 2,     delta = 1; end
    if isempty(delta), delta = 1; end
    
    sx   = sort(abs(nonzeros(x)), 'descend');
    csx  = cumsum(sx);
    nidx = find( csx - (1:numel(sx))'.*[sx(2:end); 0] >= delta ...
         + 2*eps(delta),1);
    if ~isempty(nidx)
        dx   = ( csx(nidx) - delta ) /nidx;
        xout = x.*( 1 - dx./ max(abs(x), dx) );
    else
        xout = x;
    end
    
end

% DECOPT v.1.0 by Quoc Tran-Dinh and Volkan Cevher.
% Copyright 2014 Laboratory for Information and Inference Systems (LIONS)
%                EPFL Lausanne, 1015-Lausanne, Switzerland.
% See the file LICENSE for full license information.