function [DNLmax, DNLmin, INLmax, INLmin] = Static_test(code, Nbit)
% 使用模拟法（Monte Carlo）估算正弦波概率密度，进而计算 DNL / INL

% 限制输出码范围
code = round(code);
code(code < 0) = 0;
code(code > 2^Nbit - 1) = 2^Nbit - 1;

% 获取待测 ADC 的直方图
[counts, ~] = histcounts(code, 0:2^Nbit);
M = sum(counts);  % 总采样点数

% === 模拟法计算 ideal_counts === %
% 使用高分辨率模拟一个理想正弦波输入，然后通过理想 ADC 编码再统计直方图

sim_N = 10^6;  % 模拟输入点数（越大越准）
x = linspace(0, 2*pi, sim_N);  % 相位
v = 0.5 + 0.5 * sin(x);  % 正弦输入幅值范围映射到 [0,1]

% 将输入量化为理想 ADC 输出码（无失调、理想分辨率）
ideal_code = floor(v * (2^Nbit));  % 码值范围 [0, 2^Nbit - 1]
ideal_code(ideal_code >= 2^Nbit) = 2^Nbit - 1;  % 防止上溢

% 统计理想概率密度
[ideal_counts_sim, ~] = histcounts(ideal_code, 0:2^Nbit);
ideal_pdf = ideal_counts_sim / sim_N;

% 按照实际采样点数归一化理想 bin 计数
ideal_counts = ideal_pdf * M;

% === DNL 和 INL 计算 === %
DNL = (counts - ideal_counts) ./ ideal_counts;
INL = cumsum(DNL);

% 最大/最小值输出
DNLmax = max(DNL);
DNLmin = min(DNL);
INLmax = max(INL);
INLmin = min(INL);

% === 可视化 === %
figure;
subplot(2,1,1);
bar(DNL);
xlabel('Code');
ylabel('DNL (LSB)');
title('DNL (Sinusoidal Histogram, Simulated)');
grid on;

subplot(2,1,2);
bar(INL);
xlabel('Code');
ylabel('INL (LSB)');
title('INL (Cumulative)');
grid on;

end
