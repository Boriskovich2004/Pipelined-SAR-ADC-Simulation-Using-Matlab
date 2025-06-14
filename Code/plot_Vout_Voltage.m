function plot_Vout_Voltage(Vout, Vin, Nsample, fs, L_us, R_us, Show_Input)

t = (0:Nsample-1)' / fs * 1e6; % 单位 us
index = find(t >= L_us & t <= R_us);

time = t(index);

figure;
stairs(time, Vout(index), 'b', 'LineWidth', 1.5); hold on;
if Show_Input
    plot(time, Vin(index), 'r', 'LineWidth', 1);
    legend('ADC Output', 'Input Signal');
else
    legend('ADC Output');
end
xlabel('Time (us)');
ylabel('Voltage (V)');
title('ADC Reconstructed Output Signal');
grid on;
end
