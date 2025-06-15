function [vector_Vresp, vector_Vresn] = plot_CDAC_Voltage(Vdacp, Vdacn, Vin_p, Vin_n, Nbit, Nsample, fs, L_us, R_us, Show_Input)

t = (0:Nsample*(Nbit+1)-1)' / fs * 1e6; % 单位 us
index = find(t >= L_us & t <= R_us);

Vresp = zeros(1,(Nbit+1)*Nsample);
Vresn = zeros(1,(Nbit+1)*Nsample);
for i=(1:Nsample)
    Vresp(1+(i-1)*(Nbit+1):(i)*(Nbit+1)) = Vdacp(i,:);
    Vresn(1+(i-1)*(Nbit+1):(i)*(Nbit+1)) = Vdacn(i,:);
end
% Vresp = Vdacp(:, Nbit+1);
% Vresn = Vdacn(:, Nbit+1);

vector_Vresp = Vresp(index);
vector_Vresn = Vresn(index);
time = t(index);

figure;
stairs(time, vector_Vresp, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5); hold on;
stairs(time, vector_Vresn, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
if Show_Input
    plot(time, Vin_p(index) - Vin_n(index), 'r--', 'LineWidth', 1);
    legend('CDAC Residue', 'Input Signal');
else
    legend('CDAC Residue');
end
xlabel('Time (us)');
ylabel('Voltage (V)');
title('First-stage CDAC Residue Voltage');
grid on;
end
