clc;
clear;
close all;

Nexp = 10; % 实验次数
iter_list = [50 100 150 200 250];  % 评估的不同迭代次数
avg_SNDR = zeros(size(iter_list));
avg_ENOB = zeros(size(iter_list));

set(0, 'DefaultFigureVisible', 'off');  % 关闭图像显示

for idx = 1:length(iter_list)
    iter = iter_list(idx);
    SNDR_sum = 0;
    ENOB_sum = 0;
    for k = 1:Nexp
        seed = k * 10; % 每次实验用不同种子
        [sndr_tmp, enob_tmp] = run_skew_calibration_cor(iter, seed);
        SNDR_sum = SNDR_sum + sndr_tmp;
        ENOB_sum = ENOB_sum + enob_tmp;
    end
    avg_SNDR(idx) = SNDR_sum / Nexp;
    avg_ENOB(idx) = ENOB_sum / Nexp;
end

set(0, 'DefaultFigureVisible', 'on');   % 恢复图像显示

figure;
subplot(2,1,1);
plot(iter_list, avg_SNDR, '-o', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('平均 SNDR (dB)');
title('SNDR 随迭代次数的变化');
grid on;

subplot(2,1,2);
plot(iter_list, avg_ENOB, '-s', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('平均 ENOB (bit)');
title('ENOB 随迭代次数的变化');
grid on;