function [Voutp, Voutn] = RA_FUN(Vinp, Vinn, Av, GBW, Cs, Ch, Vcm, offset, noise, fs)

Av_linear = 10^(Av/20); % dB to linear gain
beta = Cs/Ch;
Av_close = Av_linear/(1 + Av_linear*beta);
% 添加失调与热噪声
Vin_diff = Vinp - Vinn + offset + noise * randn(size(Vinp));
Vres_diff = Vin_diff * Av_close;

% 模拟放大器带宽限制影响（可选滤波建模）
tau = 1 / (2 * pi * GBW); % 时间常数
alpha = exp(-1 / (fs * tau));
Vres_filtered = filter(1 - alpha, [1, -alpha], Vres_diff);

Vres_final = Vres_filtered;

% 转为差模输出
Voutp = Vcm + Vres_final / 2;
Voutn = Vcm - Vres_final / 2;

end
