% Benchmark basisPursuitPlusL1 against decopt

addpath('../../../basisPursuitPlusL1/');
addpath('../');
addpath('../functions');
addpath('../proxs');

% =========================================================================
% Define experiments' parameters

FILENAME_RES = 'Benchm_n50_n5000_MC_50_NOCVX_mode2.mat';

% Dimensions of matrix A: m x n
n_vec = [50 100 500 1000 2000 5000];
m_vec = [30 30  50  70   150  1000];

% cardinality of x
card_x_vec = [10 10 20 30 50  200];                 

% **********************************************
% To generate prior information
% number of components in common with x
card_i_common_vec = [5 6 15 20 20 120];           

% number of components not in the support of x
card_i_rest_vec   = [5 5 5  8  30 80];           
% **********************************************

% Number of experiments ran for each dimension
MONTE_CARLO = 50;   

% If true, compares with cvx solution; otherwise, it compares with original
% (created) vector x
COMPARE_WITH_CVX = 0; 

% ADMM parameters
MAX_ITER = 5000;

% Decopt parameters
tolx        = 1e-5;
maxiter     = 5000;
verbosity   = -1;

MODE = 2;
% =========================================================================


% =========================================================================
% Experiments

% Number of experiments
num_exp = length(n_vec);

% Variables to save data
error_ADMM = zeros(num_exp, MONTE_CARLO);
error_DECO = zeros(num_exp, MONTE_CARLO);
time_ADMM  = zeros(num_exp, MONTE_CARLO);
time_DECO  = zeros(num_exp, MONTE_CARLO);

for exp_num = 1 : num_exp
    
    fprintf('Executing experiment %d (out of %d)\n', exp_num, num_exp);
    
    for mc = 1 : MONTE_CARLO
        
        n             = n_vec(exp_num);
        m             = m_vec(exp_num);
        card_x        = card_x_vec(exp_num);
        card_i_common = card_i_common_vec(exp_num);
        card_i_rest   = card_i_rest_vec(exp_num);
        
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
        b = A*x;
        % ************************************************************
        
        % ************************************************************
        % ADMM
        
        if MODE == 1
            t_ADMM_aux = cputime;
            A_pinv = pinv(A);
            [x_ADMM, k1] = basisPursuitPlusL1(b, w, 1, 1, A, A_pinv, MAX_ITER);
            t_ADMM = cputime - t_ADMM_aux;
        else
            A_handler  = @(x) A*x;
            AT_handler = @(y) A'*y;            
            t_ADMM_aux = cputime;
            [x_ADMM, k1] = basisPursuitPlusL1(b, w, 1, 2, A_handler,...
                AT_handler, MAX_ITER);
            t_ADMM = cputime - t_ADMM_aux;
        end
        % ************************************************************
        
        % ************************************************************
        % decopt 
        
        if MODE == 1
            t_DECO_aux = cputime;
            [x_DECO, output] = decoptL1L1(b, w, 1, 1, A, []);
            t_DECO = cputime - t_DECO_aux;
        else
            t_DECO_aux = cputime;
            [x_DECO, output] = decoptL1L1(b, w, 1, 2, A_handler, AT_handler);
            t_DECO = cputime - t_DECO_aux;
        end
        % ************************************************************
        
        % ************************************************************
        % Compute errors and store results
        
        % If we compare with cvx, set x_opt as its solution; otherwise,
        % set x_opt = x;
        if COMPARE_WITH_CVX
            % First execute cvx
            cvx_begin quiet
                cvx_precision high
                variable x_cvx(n);
                minimize(norm(x_cvx,1) + norm(x_cvx - w,1));
                subject to
                    A*x_cvx == b;
            cvx_end
            x_opt = x_cvx;
        else
            x_opt = x;
        end

        error_ADMM(exp_num, mc) = norm(x_ADMM - x_opt)/norm(x_opt);
        error_DECO(exp_num, mc) = norm(x_DECO - x_opt)/norm(x_opt);
        time_ADMM(exp_num, mc)  = t_ADMM;
        time_DECO(exp_num, mc)  = t_DECO;

        save(FILENAME_RES, 'error_ADMM', 'error_DECO', ...
            'time_ADMM', 'time_DECO', 'exp_num');
        % ************************************************************        
    end    
end
% =========================================================================

%%
% =========================================================================
% Print results

%load(FILENAME_RES);

error_ADMM_av = zeros(exp_num, 1);
error_DECO_av = zeros(exp_num, 1);
time_ADMM_av  = zeros(exp_num, 1);
time_DECO_av  = zeros(exp_num, 1);

for c_av = 1 : exp_num
    error_ADMM_av(c_av) = mean(error_ADMM(c_av,:)); 
    error_DECO_av(c_av) = mean(error_DECO(c_av,:)); 
    time_ADMM_av(c_av)  = mean(time_ADMM(c_av,:)); 
    time_DECO_av(c_av)  = mean(time_DECO(c_av,:)); 
end

n_vec_pl = n_vec(1:exp_num);

figure(1);clf;
semilogy(n_vec_pl, error_ADMM_av, 'bo-');
hold on;
semilogy(n_vec_pl, error_DECO_av, 'rs-');
xlabel('n');
title('Relative error (solution vector)');
legend('ADMM', 'DECO')
grid on;
drawnow;

figure(2);clf;
plot(n_vec_pl, time_ADMM_av, 'bo-');
hold on;
plot(n_vec_pl, time_DECO_av, 'rs-');
xlabel('n');
title('Time (sec)');
legend('ADMM', 'DECO')
grid on;
drawnow;
% =========================================================================
