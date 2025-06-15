function [THD, SFDR, SNR, SNDR, ENOB] = Dynamic_test(code, fs, Npoint, En_plot, wid)
% 输入向量处理
code = double(code(:)) - mean(double(code(:))); % 去直流

% 加窗处理
switch wid
    case 0, win = ones(Npoint, 1); win_name = 'None';
    case 1, win = hamming(Npoint); win_name = 'Hamming';
    case 2, win = hann(Npoint);    win_name = 'Hann';
    otherwise, error('Invalid window type');
end

y = code(1:Npoint) .* win;
Y = abs(fft(y));
Y = Y(1:Npoint/2); % 单边谱
Y_dB = 20 * log10(Y / max(Y) + eps); % dB归一化
f_axis = linspace(0, fs / 2, Npoint / 2);

% === 主频检测 ===
[peak_val, fund_idx] = max(Y(2:end)); % 跳过 DC
fund_idx = fund_idx + 1;
fund_freq = f_axis(fund_idx);

% === SFDR ===
Y_no_fund = Y;
Y_no_fund(fund_idx) = 0;
spur = max(Y_no_fund);
SFDR = 20 * log10(peak_val / spur);

% === THD（2~5次谐波）===
harmonics = 2:5;
harm_power = 0;
for h = harmonics
    harm_idx = h * fund_idx;
    if harm_idx <= Npoint / 2
        harm_power = harm_power + Y(harm_idx)^2;
        Y_no_fund(harm_idx) = 0; % 排除谐波 bin
    end
end
THD = 10 * log10(peak_val^2 / harm_power);

% === SNR ===
signal_power = peak_val^2;
noise_bins = setdiff(1:Npoint/2, [fund_idx, harmonics .* fund_idx]);
noise_power = sum(Y(noise_bins).^2);
SNR = 10 * log10(signal_power / noise_power);

% === SNDR & ENOB ===
SNDR = 10 * log10(signal_power / (harm_power + noise_power));
ENOB = (SNDR - 1.76) / 6.02;

% === 绘图 ===
if En_plot
    figure('Name', 'SAR ADC Dynamic FFT Spectrum');
    plot(f_axis / 1e6, Y_dB, 'b'); grid on; hold on;ylim([-200, 50]);
    plot(f_axis(fund_idx) / 1e6, Y_dB(fund_idx), 'ro', 'LineWidth', 1.5);
    text(f_axis(fund_idx)/1e6, Y_dB(fund_idx) + 3, ...
        sprintf('Peak @ %.2f MHz', fund_freq / 1e6), 'FontSize', 10);

    xlabel('Frequency (MHz)');
    ylabel('Magnitude (dB)');
    title(['FFT Spectrum (', num2str(Npoint), 'pt, ', win_name, ' Window)']);

    % 显示动态指标
    txt = {
        sprintf('SNR   = %.2f dB', SNR)
        sprintf('SNDR  = %.2f dB', SNDR)
        sprintf('SFDR  = %.2f dB', SFDR)
        sprintf('THD   = %.2f dB', THD)
        sprintf('ENOB  = %.2f bit', ENOB)
        };
    text(0.7 * max(f_axis)/1e6, max(Y_dB)-5, txt, ...
        'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
end
end
