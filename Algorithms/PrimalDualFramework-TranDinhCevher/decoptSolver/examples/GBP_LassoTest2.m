% TEST OUR SOLVER FOR THE BASIS PURSUIT PROBLEM OF THE FORM:
%
%               min |w.*x|_1  s.t. A*x = b.
%
% Date: 31.12.2013
%

plist         = [5:20, 5:20,21:25];

for pid = 16:length(plist)

scale   = plist(pid);
n       = scale*1024;
m       = scale*341; 
ng      = round(m/4);
    
% Load the data ...
load_file = ['p_5noise_5%_', num2str(m), '_', num2str(n), '_', num2str(ng), '_', num2str(pid), '.mat'];
load(load_file);
    
%% Call the solver (Algorithm 2).
param.Algorithm     = 2;
param.adaptStepSize = 0;
time1b       = tic;
[x1b, out1b] = decoptSolver('GBP', A, b, param, 'x0', x0, 'groups', groups);
time1b       = toc(time1b);

% Evaluate the objective values and feasibility gap.
fx1b   = sum(weights.*sqrt(mgroups*x1b.^2));
feas1b = norm(A*x1b - b, 2); 

% %% If exist YALL1, then compare with it.
% if isExistSGPL1 
%     tolx2 = 1e-8;
%     fprintf('---> SPGL1 ...\n');
%     opts4 = spgSetParms('verbosity',  2, ... %isprint, ...
%                         'weights',    weights, ...
%                         'bpTol',      1e-2*tolx2, ...
%                         'lsTol',      1e-2*tolx2, ...
%                         'optTol',     tolx2, ...
%                         'decTol',     tolx2, ...
%                         'iterations', maxiters);
%     sigma_b = sigma*norm(b);
%     time4             = tic;
%     [x4, R, G, info4] = spg_group(A, b, gindices, sigma_b, opts4);
%     time4             = toc(time4);
%     
%     % Evaluate the objective values and feasibility gap.
%     fx4   = sum(weights.*sqrt(mgroups*x4.^2));
%     feas4 = norm(A*x4 - b, 2); 
% end

%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM(1): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1a, feas1a/norm(b));
fprintf('+ DECOM(1): Iterations: %4d, Time(s) = %3.4f\n', out1a.iter, out1a.total_time);
fprintf('+ DECOM(1): Number of Ax and ATy are %4d and %4d\n', out1a.cntA, out1a.cntAt);
fprintf('+ DECOM(1): Reconvery error: %4.12f\n', norm(x1a - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(2): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1b, feas1b/norm(b));
fprintf('+ DECOM(2): Iterations: %4d, Time(s) = %3.4f\n', out1b.iter, out1b.total_time);
fprintf('+ DECOM(2): Number of Ax and ATy are %4d and %4d\n', out1b.cntA, out1b.cntAt);
fprintf('+ DECOM(2): Reconvery error: %4.12f\n', norm(x1b - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(3): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1c, feas1c/norm(b));
fprintf('+ DECOM(3): Iterations: %4d, Time(s) = %3.4f\n', out1c.iter, out1c.total_time);
fprintf('+ DECOM(3): Number of Ax and ATy are %4d and %4d\n', out1c.cntA, out1c.cntAt);
fprintf('+ DECOM(3): Reconvery error: %4.12f\n', norm(x1c - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(4): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1d, feas1d/norm(b));
fprintf('+ DECOM(4): Iterations: %4d, Time(s) = %3.4f\n', out1d.iter, out1d.total_time);
fprintf('+ DECOM(4): Number of Ax and ATy are %4d and %4d\n', out1d.cntA, out1d.cntAt);
fprintf('+ DECOM(4): Reconvery error: %4.12f\n', norm(x1d - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(5): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1e, feas1e/norm(b));
fprintf('+ DECOM(5): Iterations: %4d, Time(s) = %3.4f\n', out1e.iter, out1e.total_time);
fprintf('+ DECOM(5): Number of Ax and ATy are %4d and %4d\n', out1e.cntA, out1e.cntAt);
fprintf('+ DECOM(5): Reconvery error: %4.12f\n', norm(x1e - x_org)/max(norm(x_org), 1));

fprintf('+ DECOM(6): BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx1f, feas1f/norm(b));
fprintf('+ DECOM(6): Iterations: %4d, Time(s) = %3.4f\n', out1f.iter, out1f.total_time);
fprintf('+ DECOM(6): Number of Ax and ATy are %4d and %4d\n', out1f.cntA, out1f.cntAt);
fprintf('+ DECOM(6): Reconvery error: %4.12f\n', norm(x1f - x_org)/max(norm(x_org), 1));

if isExistYALL1 
    fprintf('+ YALL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx3, feas3/norm(b));
    fprintf('+ YALL1: Iterations: %4d, Time(s) = %3.4f\n', out3.iter, time3);
%    fprintf('+ YALL1: Number of Ax and ATy are %4d and %4d\n', out3.cntA, out3.cntAt);
    fprintf('+ YALL1: Reconvery error: %4.12f\n', norm(x3 - x_org)/max(norm(x_org), 1));
end
if isExistSGPL1 
    fprintf('+ SPGL1: BP-problem: [f(x), |A*x-b|/|b|] = [%3.12f, %5.12f]\n', fx4, feas4/norm(b));
    fprintf('+ SPGL1: Iterations: %4d, Time(s) = %3.4f\n', info4.iter, info4.timeTotal);
    fprintf('+ SPGL1: Number of Ax and ATy are %4d and %4d\n', info4.nProdA, info4.nProdAt);
    fprintf('+ SPGL1: Reconvery error: %4.12f\n', norm(x4 - x_org)/max(norm(x_org), 1));
end

%%
clear A;
save(['modified/', save_file], '-v7.3');

end
%% END OF THE TEST.