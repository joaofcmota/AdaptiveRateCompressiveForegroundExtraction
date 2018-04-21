% PURPOSE: This function tests the logistic regression solver.
%
% DATE: 08.04.2013.
%
% INPUTS:   file_name  : the name of the data file.
%           method     : proximal-gradient or proximal-Newton method.
%

% Data ...
files_train     = {'a2a.txt'};
files_test      = {'a2a.t'};
featuresNo      = [123];
number_of_tests = length(featuresNo);


% Set of regularization parameters.
sC              = 10.^(-3:0.5:3);
number_of_Cs    = length(sC);

% Composite function optimization parameters
RelTolX         = 1e-9; % Stopping criteria tolerance

% Kernel bias
bias            = 0;

% Initialization of statistics variables
number_tests_t  = number_of_Cs; % Number of test for time for values in the subset sC
% For composite function minization
stat_acc        = zeros(number_of_tests,number_of_Cs);
stat_obj        = zeros(number_of_tests,number_of_Cs);
stat_time       = cell(number_of_tests,number_tests_t);

% For LIBSVM
stat_acc_lib        = zeros(number_of_tests,number_of_Cs);
stat_obj_lib        = zeros(number_of_tests,number_of_Cs);
stat_time_lib       = cell(number_of_tests,number_tests_t);

% To plot time and accuracy for LIBSVM we change the tolerance, making the 
% algorithm terminate early and measure the objective and time
% We fix a high tolerance to time memeory transfer operations and other
% non related to the algorithm processes
tol_time            = linspace(1,10^-5,20);

% tol_time_lib        = [10^3 tol_time];
tol_time_lib = [10 7.5 5  2 1 0.75 0.5 0.25 linspace(10^-1,10^-5,5)];

Cs_time = zeros(number_of_tests,number_tests_t);

for tt=1:number_of_tests
    
  file_name       = files{tt};
  if ~isempty(files_test{tt})
      file_name_test = files_test{tt};
  end
  p               = featuresNo(tt);
% file_name       = 'a3a';
% file_name_test  = 'a3a_t'; % Comment this to switch to single file mode

% file_name       = ['a',num2str(tt),'a'];
% file_name_test  = ['a',num2str(tt),'a_t'];
% p               = 123;     % Number of features - defined for every set
%                            % Trailing zeros are ignored in the sparse
%                            % description of the input
% 
separation      = 0.7;     % Percent of data atributed to training - SFM
 

SFM             = ~exist('file_name_test','var');

if SFM == 1
    %-----------------------------
    % Single file mode
    %-----------------------------
    
    %data_path = which('logistic_exam01.m');
    data_path = '/Users/quoctd/Data/MyPostDoc/Epfl-SVN2013/matlab/SCOPT-v1.0/examples/logistic_data/';
    %dirData   = [data_path(1:end-17), 'logistic_data/'];
    dirData   = data_path;

    % Read data from the data set.
    if exist('libsvmread.c','file')
    [yhat, What] = libsvmread(fullfile(dirData, file_name));
    %save([dirData, mat_file_name], 'yhat', 'What');
    else
        mdata = load([dirData, mat_file_name]);
        yhat  = mdata.yhat;
        What  = mdata.What;
    end


    if isempty(yhat) || isempty(What)
        error('Can not read the data from the file!');
    end
    
    % Split data into training and test - not randomized for now
    N = size(What,1);

    size_training = round(0.7*N);
    size_test     = N - size_training;

    Wtest = What(size_training+1:end,:);
    ytest = yhat(size_training+1:end,:);
    
    What = What(1:size_training,:);
    yhat = yhat(1:size_training,:);
    
     fprintf(' + We are solving problem: %s (p = %d, N_training = %d, N_test = %d)\n', ...
         file_name, p, size_training,size_test);
    
else
    %-----------------------------
    % Multi file mode
    %-----------------------------
    

    % Training set
    %data_path = which('logistic_exam01.m');
    %dirData   = [data_path(1:end-17), 'logistic_data/'];
    data_path = '/Users/quoctd/Data/MyPostDoc/Epfl-SVN2013/matlab/SCOPT-v1.0/examples/logistic_data/';
    dirData   = data_path;

    % Read data from the data set.
    if exist('libsvmread.c','file')
    [yhat, What] = libsvmread(fullfile(dirData, file_name));
    %save([dirData, mat_file_name], 'yhat', 'What');
    else
        mdata = load([dirData, mat_file_name]);
        yhat  = mdata.yhat;
        What  = mdata.What;
    end


    if isempty(yhat) || isempty(What)
        error('Can not read the data from the file!');
    end
    
    
    
    % Test set
    %data_path = which('logistic_exam01.m');
    %dirData   = [data_path(1:end-17), 'logistic_data/'];
    data_path = '/Users/quoctd/Data/MyPostDoc/Epfl-SVN2013/matlab/SCOPT-v1.0/examples/logistic_data/';
    dirData   = data_path;

    % Read data from the data set.
    if exist('libsvmread.c','file')
    [ytest, Wtest] = libsvmread(fullfile(dirData, file_name_test));
    %save([dirData, mat_file_name], 'yhat', 'What');
    else
        
    end


    if isempty(ytest) || isempty(Wtest)
        error('Can not read the data from the file!');
    end
    
    % Padding with zeros
    size_training   = size(What,1);
    size_test       = size(Wtest,1);
    
    p1              = size(What,2);
    p2              = size(Wtest,2);

    if p1<p
        %Padd the training set with zeros
        What = [What zeros(size_training,p-p1)];
    end
    
    if p2<p
        %Padd the test set with zeros
        Wtest = [Wtest zeros(size_test,p-p2)];
    end
    
    fprintf(' + We are solving problem: training %s - test %s (p = %d, N_training = %d, N_test = %d)\n', ...
         file_name,file_name_test, p, size_training,size_test);
end




%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
% What - feature matrix for training
% yhat - decision for training
% Wtest - feature matrix for testing
% ytest - decision for testing
% size_training - number of training data
% size_test - number of test data
% p - dimension of features

 

% First tests are to explore relation to parameter C


for cc=1:length(sC)
    
C = sC(cc);


%% Call the decomposition algorithm.
fprintf('------------------------------------------\n');
fprintf('     Decomposition algorithm              \n');
fprintf('------------------------------------------\n');

% Parameters
opts0 = decompParamSettings([]);
opts0.MaxIters = 10000;
opts0.RelTolX  = RelTolX;
opts0.isNormL1 = useL1;
opts0.isFixedBias = false;



% Generate the initial point.
x0    = 0*ones(p, 1);



% To match parameters in LIBSVM we need 2*C as the regularization param for
% L2 norm

%% Set the parameters.
param.MaxIters      = 10000;
param.Verbosity     = 2;
param.RelTolX       = RelTolX*1e-2;
param.RelTolFeas    = RelTolX*1e-2;
param.saveHistMode  = 0;
param.Algorithm     = 1;
param.InnerMaxIters = 50;
param.adaptStepSize = 0;

%% User-define prox-functions.
rho = 1/C;
if 0%isL1Norm
    proxOpers{1} = @(x, gamma)([proxL1norm(x(1:end-1), gamma); x(end)]);
    proxOpers{3} = @(x)( norm(rho.*x(1:end-1), 1) );
else
    proxOpers{1} = @(x, gamma)([proxSquareL2norm(x(1:end-1), gamma); x(end)]);
    proxOpers{3} = @(x)(0.5*sum(rho.*x(1:end-1).^2));
end
proxOpers{2} = @(x, gamma)(proxHingeLoss(x, gamma, yhat, yhat.^2));
proxOpers{4} = @(x)(sum(max(1 - yhat.*x, 0)));

% Regenerate the data.
[N, p] = size(What);
What2 = [What, ones(N,1)];
bias2 = zeros(size(bias));
x02   = [x0; 0];
lbx   = [-inf*ones(p,1); -10];
ubx   = [+inf*ones(p,1); +10];

% Call the solver.
[xsol_dc2f, output_dc] = decoptSolver('UserDef', What2, bias2, param, ...
                                      'x0', x02, 'RegPar', rho, ...
                                      'Label', yhat, 'Prox', proxOpers, ...
                                      'lbx', lbx, 'ubx', ubx);    
xsol_dc2 = xsol_dc2f(1:end-1);
new_bias = xsol_dc2f(end);
fsol.x   = xsol_dc2;
fsol.bias = new_bias;
output_dc.time = output_dc.total_time;

% Model training
%[fsol, output_dc]           = hingelossMainAlg2(What, yhat, 0, ...
%                             1/C, opts0);

% Testing

[fval_dc,fxval,gxval,acc]   = scopt_hinge_fx_eval(What, fsol.x, fsol.bias,... 
                              C, yhat,useL1);

[fval_test,fxval,gxval,acc_test]   = scopt_hinge_fx_eval(Wtest, fsol.x, fsol.bias,... 
                             C, ytest,useL1);

nnz_dc     = nnz(round(fsol.x*1e10)*1e-10);
fprintf('**** FINAL OUTPUTS OF DECOMPOSITION ALGORITHM ****\n');
fprintf(' + Number of iterations: %5d\n', output_dc.iter);
fprintf(' + The objective value for test data: %6.15f\n', fval_dc);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*fxval,0.5*gxval);
fprintf(' + Accuracy %1.5f\n',acc);
fprintf(' + CPU time: %6.4f\n', output_dc.time);
fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_dc);
fprintf('********************************************************\n');

stat_acc(tt,cc)        = acc_test;
stat_obj(tt,cc)        = fval_dc;

%%
% %---------------------------------------------------
% % Compare with LIBSVM
% %---------------------------------------------------


if ~useL1
    
    
% %Here we determine the parameters
%  svmstruct.Bias  = bias;
%  opt             = statset('Display','iter','MaxIter',15000);
% model           = svmtrain(yhat,What ,...
%                     'kernel_function','linear',...
%                   'method','SMO','options',opt,'tolkkt',10^-3,...
%                   'kktviolationlevel',0,'boxconstraint',C,...
%                   'autoscale',true,'SVMSTRUCT',svmstruct);


%Description of paramenteres for LIBSVM
% libsvm_options:
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC		(multi-class classification)
% 	1 -- nu-SVC		(multi-class classification)
% 	2 -- one-class SVM
% 	3 -- epsilon-SVR	(regression)
% 	4 -- nu-SVR		(regression)
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_instance_matrix)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n : n-fold cross validation mode
% -q : quiet mode (no outputs)



%Cross-validation



options_libsvm = ['-s 0 -t 0 -e 0.001 -c ',num2str(C)];
    

 tic;
 model           = svmtrain(yhat,What ,options_libsvm);
 time_test = toc;
% % Here we compute the accuracy              
[predicted_label, accuracy, decision_values] = ...
     svmpredict(ytest, Wtest, model);
 
 % Parameters from the model
 W_libsvm = model.SVs' * model.sv_coef;
 bias_libsvm = -model.rho;
 
 [fval_svm, fxval, gxval, acc]   = scopt_hinge_fx_eval(What, W_libsvm,  bias_libsvm,... 
                             C, yhat,useL1);
 
 [fval_test, fxval, gxval, acc_test]   = scopt_hinge_fx_eval(Wtest, W_libsvm,  bias_libsvm,... 
                             C, ytest,useL1);
                         
fprintf(' + The objective value for test data: %6.15f\n', fval_svm);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*fxval,0.5*gxval);
fprintf(' + Accuracy %1.5f\n',acc);

stat_acc_lib(tt,cc)        = acc_test;
stat_obj_lib(tt,cc)        = fval_svm;
end





end

% This tests are to measure time

a = sortrows([[(stat_obj(tt,:)-stat_obj_lib(tt,:))./stat_obj(tt,:)]',[1:number_of_Cs]']);
Cs_time(tt,:) = sC(a(1:number_tests_t,2));

for cc=1:number_tests_t
    C = Cs_time(tt,cc);
    
    
%% Call the decomposition algorithm.
fprintf('------------------------------------------\n');
fprintf('     Decomposition algorithm              \n');
fprintf('------------------------------------------\n');

% Parameters
opts0 = decompParamSettings([]);
opts0.MaxIters = 10000;
opts0.RelTolX  = RelTolX;
opts0.isNormL1 = useL1;
opts0.isFixedBias = false;
opts0.saveHistMode = 2;
opts0.Verbosity = 0;


% Generate the initial point.
x0    = 0*ones(p, 1);



% To match parameters in LIBSVM we need 2*C as the regularization param for
% L2 norm

% Model training
[fsol, output_dc]           = hingelossMainAlg2(What, yhat, 0, ...
                             1/C, opts0);

% Testing

[fval_dc,fxval,gxval,acc]   = scopt_hinge_fx_eval(What, fsol.x, fsol.bias,... 
                              C, yhat,useL1);



nnz_dc     = nnz(round(fsol.x*1e10)*1e-10);
fprintf('**** FINAL OUTPUTS OF DECOMPOSITION ALGORITHM ****\n');
fprintf(' + Number of iterations: %5d\n', output_dc.iter);
fprintf(' + The objective value for test data: %6.15f\n', fval_dc);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*fxval,0.5*gxval);
fprintf(' + Accuracy %1.5f\n',acc);
fprintf(' + CPU time: %6.4f\n', output_dc.time);
fprintf(' + Number of nonzero elements of xsol: %d\n', nnz_dc);
fprintf('********************************************************\n');

stat_time{tt,cc} = [output_dc.hist.time, output_dc.hist.fx_val];


time_vs_obj_lib = zeros(length(tol_time_lib),2);
    
for kk=1:length(tol_time_lib);
 
tolerance_lib = tol_time_lib(kk);    
%%
% %---------------------------------------------------
% % Compare with LIBSVM
% %---------------------------------------------------


if ~useL1
    
    
% %Here we determine the parameters
%  svmstruct.Bias  = bias;
%  opt             = statset('Display','iter','MaxIter',15000);
% model           = svmtrain(yhat,What ,...
%                     'kernel_function','linear',...
%                   'method','SMO','options',opt,'tolkkt',10^-3,...
%                   'kktviolationlevel',0,'boxconstraint',C,...
%                   'autoscale',true,'SVMSTRUCT',svmstruct);


%Description of paramenteres for LIBSVM
% libsvm_options:
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC		(multi-class classification)
% 	1 -- nu-SVC		(multi-class classification)
% 	2 -- one-class SVM
% 	3 -- epsilon-SVR	(regression)
% 	4 -- nu-SVR		(regression)
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_instance_matrix)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n : n-fold cross validation mode
% -q : quiet mode (no outputs)



%Cross-validation



options_libsvm = ['-s 0 -t 0 -e ',num2str(tolerance_lib,'%f'),' -c ',num2str(C)];
    

 tic;
 model           = svmtrain(yhat,What ,options_libsvm);
 time_test = toc;
% % Here we compute the accuracy              
[predicted_label, accuracy, decision_values] = ...
     svmpredict(ytest, Wtest, model);
 
 % Parameters from the model
 W_libsvm = model.SVs' * model.sv_coef;
 bias_libsvm = -model.rho;
 
 [fval_svm, fxval, gxval, acc]   = scopt_hinge_fx_eval(What, W_libsvm,  bias_libsvm,... 
                             C, yhat,useL1);
 
                         
fprintf(' + The objective value for test data: %6.15f\n', fval_svm);
fprintf(' + As %6.15f and %6.15f (for penalty)\n',C*fxval,0.5*gxval);
fprintf(' + Accuracy %1.5f\n',acc);

time_vs_obj_lib(kk,1) = time_test;
time_vs_obj_lib(kk,2) = fval_svm;
end    
    
end

% time_vs_obj_lib(:,1) = time_vs_obj_lib(:,1) - time_vs_obj_lib(1,1);  
 
stat_time_lib{tt,cc}=time_vs_obj_lib;    
    
    
    
    
end



%Generating plots
figure();
semilogx(sC,stat_acc(tt,:),'r-');
xlabel('C');
ylabel('Accuracy');
hold on
semilogx(sC,stat_acc_lib(tt,:),'b--');
hold off;

minimum_value = min([stat_obj(tt,:)]);
figure();
semilogx(sC,(stat_obj(tt,:))./stat_obj_lib(tt,:),'r-');
xlabel('C');
ylabel('Objective value');
hold on
semilogx(sC,ones(size(sC)),'b--');
hold off;
% semilogx(sC,stat_obj_lib(tt,:)/minimum_value,'bx--');
% hold off;

% This does not work
% for cc=1:number_tests_t
% figure();
% semilogy(stat_time{cc}(:,cc),stat_time{cc}(:,2),'r-');
% hold on;
% semilogy(stat_time_lib{cc}(:,1),stat_time_lib{cc}(:,2),'bx');
% hold off;
% end


spath = which('logistic_exam02.m');
spath = spath(1:end-length('logistic_exam02.m'));
save([spath, 'Figs/prob_',file_name,'_C', num2str(min(sC)),'_',num2str(max(sC)),'_Tol',num2str(RelTolX)],'stat_acc','stat_acc_lib','stat_obj','stat_obj_lib','sC','RelTolX');

end





