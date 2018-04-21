clear all; clc;

C = 1000;
DATA_FILE = 'a1a.txt';
nprob = 10;

%DATA_FILE = 'news20.binary_new';
%nprob  = 6;

for pk = 1:nprob
    
    file_names = [DATA_FILE, '_', num2str(pk), '_', num2str(C), '.mat'];
    DATA = load(file_names);
    
    % Objective.
    fx_trn_dc(pk, :) = DATA.fx_list(1, :);
    fx_tst_dc(pk, :) = DATA.fx_list(2, :);
    fx_trn_sv(pk, :) = DATA.fx_list(3, :);
    fx_tst_sv(pk, :) = DATA.fx_list(4, :);
    
    % Accuracy.
    ac_trn_dc(pk, :) = DATA.acc_list(1, :);
    ac_tst_dc(pk, :) = DATA.acc_list(2, :);
    ac_trn_sv(pk, :) = DATA.acc_list(3, :);
    ac_tst_sv(pk, :) = DATA.acc_list(4, :);
    
    % Time
    tm_trn_dc(pk, :) = DATA.time_list(1, :);
    tm_trn_sv(pk, :) = DATA.time_list(2, :);
end

%%
sC                   = DATA.sC;
nrn                  = nprob;
fx_list_avg_dc       = sum(fx_trn_dc)/nrn;
fx_list_min_dc       = abs( min(fx_trn_dc) - fx_list_avg_dc );
fx_list_max_dc       = abs(max(fx_trn_dc) - fx_list_avg_dc);

fx_list_avg_sv       = sum(fx_trn_sv)/nrn;
fx_list_min_sv       = abs( min(fx_trn_sv) - fx_list_avg_sv);
fx_list_max_sv       = abs( max(fx_trn_sv) - fx_list_avg_sv);

fx_list_avg_dc_test  = sum(fx_tst_dc)/nrn;
fx_list_min_dc_test  = abs( min(fx_tst_dc) - fx_list_avg_dc_test);
fx_list_max_dc_test  = abs( max(fx_tst_dc) - fx_list_avg_dc_test);

fx_list_avg_sv_test  = sum(fx_tst_sv)/nrn;
fx_list_min_sv_test  = abs( min(fx_tst_sv) - fx_list_avg_sv_test);
fx_list_max_sv_test  = abs( max(fx_tst_sv) - fx_list_avg_sv_test);

acc_list_avg_dc      = sum(ac_trn_dc)/nrn;
acc_list_min_dc      = abs( min(ac_trn_dc) - acc_list_avg_dc);
acc_list_max_dc      = abs( max(ac_trn_dc) - acc_list_avg_dc);

acc_list_avg_sv      = sum(ac_trn_sv)/nrn;
acc_list_min_sv      = abs( min(ac_trn_sv) - acc_list_avg_sv);
acc_list_max_sv      = abs( max(ac_trn_sv) - acc_list_avg_sv);

acc_list_avg_dc_test = sum(ac_tst_dc)/nrn;
acc_list_min_dc_test = abs( min(ac_tst_dc) - acc_list_avg_dc_test);
acc_list_max_dc_test = abs( max(ac_tst_dc) - acc_list_avg_dc_test);

acc_list_avg_sv_test = sum(ac_tst_sv)/nrn;
acc_list_min_sv_test = abs( min(ac_tst_sv) - acc_list_avg_sv_test);
acc_list_max_sv_test = abs( max(ac_tst_sv) - acc_list_avg_sv_test);

time_list_avg_dc     = sum(tm_trn_dc)/nrn;
time_list_min_dc     = abs( min(tm_trn_dc) - time_list_avg_dc);
time_list_max_dc     = abs( max(tm_trn_dc) - time_list_avg_dc);

time_list_avg_sv     = sum(tm_trn_sv)/nrn;
time_list_min_sv     = abs( min(tm_trn_sv) - time_list_avg_sv);
time_list_max_sv     = abs( max(tm_trn_sv) - time_list_avg_sv);

%%
figure('Name', 'Fx for training data');
subplot(3, 2, 1);
%plot(fx_list_avg_dc, 'r--'); hold on;
%plot(fx_list_avg_sv, 'b:');  
errorbar(sC, fx_list_avg_dc, fx_list_min_dc, fx_list_max_dc, '-bs'); hold on;
errorbar(sC, fx_list_avg_sv, fx_list_min_sv, fx_list_max_sv, '-ro'); hold on;
title('The objective values for the training data');
legend('1P2D', 'LIBSVM');

subplot(3, 2, 2);
%plot(fx_list_avg_dc_test, 'r--'); hold on;
%plot(fx_list_avg_sv_test, 'b:');  
errorbar(sC, fx_list_avg_dc_test, fx_list_min_dc_test, fx_list_max_dc_test, '-bs'); hold on;
errorbar(sC, fx_list_avg_sv_test, fx_list_min_sv_test, fx_list_max_sv_test, '-ro'); hold on;
title('The objective values for the test data');
legend('1P2D', 'LIBSVM');

subplot(3, 2, 3);
%plot(acc_list_avg_dc, 'r--'); hold on;
%plot(acc_list_avg_sv, 'b:');  
errorbar(sC, acc_list_avg_dc, acc_list_min_dc, acc_list_max_dc, '-bs'); hold on;
errorbar(sC, acc_list_avg_sv, acc_list_min_sv, acc_list_max_sv, '-ro'); hold on;
title('The accuracy values for the training data');
legend('1P2D', 'LIBSVM');

subplot(3, 2, 4);
%plot(acc_list_avg_dc_test, 'r--'); hold on;
%plot(acc_list_avg_sv_test, 'b:');  
errorbar(sC, acc_list_avg_dc_test, acc_list_min_dc_test, acc_list_max_dc_test, '-bs'); hold on;
errorbar(sC, acc_list_avg_sv_test, acc_list_min_sv_test, acc_list_max_sv_test, '-ro'); hold on;
title('The accuracy values for the test data');
legend('1P2D', 'LIBSVM');

subplot(3, 2, 5);
%plot(time_list_avg_dc, 'r--'); hold on;
%plot(time_list_avg_sv, 'b:');  
errorbar(sC, time_list_avg_dc, time_list_min_dc, time_list_max_dc, '-bs'); hold on;
errorbar(sC, time_list_avg_sv, time_list_min_sv, time_list_max_sv, '-ro'); hold on;
title('The cpu time');
legend('1P2D', 'LIBSVM');