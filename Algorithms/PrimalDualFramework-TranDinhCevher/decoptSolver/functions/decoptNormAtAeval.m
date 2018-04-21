%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: l2nrm = decoptNormAtAeval(mode, A, maxiters, tol, AT, nx)
%%% PURPOSE: Compute the l2-norm of A^T*A by using Power method.
%%%
%%% INFORMATION:
%%%    By Quoc Tran Dinh, Laboratory for Informations and Inference Systems,
%%%       (LIONS), EPFL, Lausanne, Switzerland.
%%%    Joint work with Volkan Cevher and Anastasios Kyrillidis.
%%%    Date: 03.06.2013.
%%%    Last modified: 20.12.2013.
%%%    Contact: quoc.trandinh@epfl.ch
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l2nrm = decoptNormAtAeval(mode, A, maxiters, tol, AT, nx)

% For numeric case.
if strcmpi(mode, 'numeric')
    l2nrm = powerIterForMat(A, maxiters, tol);
    
% For operation case.
else
    if nargin < 6, error('At least six inputs are required!'); end
    l2nrm = powerIterForOper(A, AT, nx, maxiters, tol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: l2nrm = powerIterForMat(A, maxiters, tolx)
%%% PURPOSE:  Compute the norm of matrix (A'*A).
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l2nrm = powerIterForMat(A, maxiters, tol)

% Define the conjugate operator.
nx = size(A, 2);

% Initialization.
z0     = ones(nx, 1);
q      = z0/norm(z0, 2); 
relres = tol + 1; 
z      = A'*(A*q);

% The main loop.
for iter = 1:maxiters
    if relres <= tol, break; end
    
    % Compute the new iteration.
    q       = z/norm(z, 2); 
    z       = A'*(A*q);
    lambda  = q'*z; 
    z2      = conj(z);
    q2      = z2/norm(z2); 
    y1      = q2; 
    cosqy   = abs(y1'*q);
    
    if cosqy >= 5e-2
        relres = norm(z - lambda*q, 2)/cosqy;
    end
end

% Compute the norm of A'*A.
l2nrm = abs(lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: l2nrm = powerIterForOper(A, AT, nx, maxiters, tolx)
%%% PURPOSE:  Compute the norm of the operator (AT o A).
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l2nrm = powerIterForOper(A, AT, nx, maxiters, tol)

% Initialization.
z0     = ones(nx, 1);
q      = z0/norm(z0, 2); 
relres = tol + 1; 
z      = AT(A(q));

% The main loop.
for iter = 1:maxiters
    if relres <= tol, break; end
    
    % Compute the new iteration.
    q       = z/norm(z, 2); 
    z       = AT(A(q));
    lambda  = q'*z; 
    z2      = conj(z);
    q2      = z2/norm(z2); 
    y1      = q2; 
    cosqy   = abs(y1'*q);
    
    if cosqy >= 5e-2
        relres = norm(z - lambda*q, 2)/cosqy;
    end
end

% Compute the norm of A'*A.
l2nrm = abs(lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF THE IMPLEMENTATION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%