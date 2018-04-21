% Tests the decopt solver in the problem
%
%   minimize    ||x||_1 + ||x - w||_1
%      x
%   subject to  ||A*x - b|| <= sigma

addpath('../');
addpath('../functions');
addpath('../proxs');

% =========================================================================
% Define experiments' parameters

% Dimensions of matrix A: m x n
n = 16384;
m = 500;

% cardinality of x
card_x = 20;

% **********************************************
% To generate prior information
% number of components in common with x
card_i_common = 10;

% number of components not in the support of x
card_i_rest   = 5;
% **********************************************

sigma = 0.1;

% Decopt parameters
tolx        = 1e-5;
maxiter     = 5000;
verbosity   = 0;

MODE = 1;
% =========================================================================


% =========================================================================
% Experiments


% ************************************************************
% Generate data

% Generate x
x_aux = [randn(card_x,1); zeros(n-card_x,1)];
permutation_x = randperm(n);
x = x_aux(permutation_x);

% Generate w
i_aux = [randn(card_i_rest,1); zeros(n-card_i_rest,1)];
permutation_rest = randperm(n);
i = i_aux(permutation_rest);
vec_aux = [randn(card_i_common,1); zeros(n-card_i_common,1)];
vec_perm = vec_aux(permutation_x);
i = i + vec_perm;
w = x + i;

% Generate measurement matrix and vector of measurements
A = randn(m,n);
b = A*x + (sigma/m)*randn(m,1);
% ************************************************************
%%
% ************************************************************
% decopt
tic
if MODE == 1
    t_DECO_aux = cputime;
    [x_DECO, output] = decoptL1L1Noise(b, sigma, w, 1, 1, A, [], maxiter);
    t_DECO = cputime - t_DECO_aux;
else
    t_DECO_aux = cputime;
    [x_DECO, output] = decoptL1L1Noise(b, sigma, w, 1, 2, A_handler, AT_handler, maxiter);
    t_DECO = cputime - t_DECO_aux;
end
toc
% ************************************************************
%%
% % ************************************************************
% % Compare with CVX
% 
% % First execute cvx
% cvx_begin
%     cvx_precision high
%     variable x_cvx(n);
%     minimize(norm(x_cvx,1) + norm(x_cvx - w,1));
%     subject to
%         norm(A*x_cvx - b) <= sigma;
% cvx_end
% x_opt = x_cvx;
% % =========================================================================
% 
% fprintf('||x_DECO - x_opt||/||x_opt||    = %f\n', norm(x_DECO - x_opt)/norm(x_opt));
fprintf('||A*x_DECO - b||                = %f\n', norm(A*x_DECO - b));
% fprintf('||A*x_opt  - b||                = %f\n', norm(A*x_opt - b));
fprintf('||x_DECO||_1 + ||x_DECO - w||_1 = %f\n', norm(x_DECO,1) + norm(x_DECO - w,1));
% fprintf('||x_opt||_1  + ||x_opt  - w||_1 = %f\n', norm(x_opt ,1) + norm(x_opt  - w,1));