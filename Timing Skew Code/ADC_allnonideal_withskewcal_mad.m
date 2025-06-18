clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ADC Parameter %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 7;          % stage1 SAR
N2 = 9;          % stage2 SAR
N = N1 + N2 - 1; % resolution (one bit redundancy)
fs = 100e6; % Sample Rate(MHz)
ts = 1/(fs);
Vref = 1.8; % Reference Voltage
Vcm = Vref/2; % Common Mode Voltage
LSB = Vref/(2^N);
n_ch = 4; % number of channel

%%%%%%%%%% Mismatch & noise %%%%%%%%%%%
% Inter-Channel Mismatch
Mis_OS = 0.5 * LSB * randn(1,n_ch);       % Inter-channel Offset
Mis_Gain = 1 + 0.2 * LSB * randn(1,n_ch); % Inter-channel Gain Mismatch
rng(10);
Mis_TS = 0.05 * 1/fs *  randn(1, n_ch);    % Inter-channel Clock Skew

% SAR Capacitor Mismatch
sigmaC = 0.00;    % Global Capacitor Mismatch
Cp_p1  = 5e-15; % Stage1 Positive Plate Parasitic Capacitor (F)
Cp_n1  = 5e-15; % Stage1 Negative Plate Parasitic Capacitor (F)
Cp_p2  = 5e-15; % Stage2 Positive Plate Parasitic Capacitor (F)
Cp_n2  = 5e-15; % Stage2 Negative Plate Parasitic Capacitor (F)

% ResAmp Mismatch & noise
Av         = 80;    % Open-loop Gain (dB)
GBW        = 0.1e9;   % Gain-Bandwidth Product (GHz)
Amp_offset = 1e-3; % Amplifier Offset (V)
Amp_noise  = 000e-3; % Amplifier RMS Voltage Noise (V)

% SAR Comparator Mismatch & noise
Comp_offset1 = 5e-3;% Stage1 SAR Comp Offset (V)
Comp_offset2 = 5e-3;% Stage2 SAR Comp Offset (V)
Comp_noise   = 000e-3;% SAR Comp RMS Voltage Noise (V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Environment Parameter %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.38e-23;       % Bolzman Constant
T = 300;            % Temperature (K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sub ADC Parameter %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_sub = fs / n_ch;    % Sub ADC Sample Rate
ts_sub = 1 / (fs_sub); % Sub ADC Sample Cycle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% input signal %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 2^22;            % Sample Points
Nsample = num;         % FFT Points
Vfs = Vref;            % Input Signal Full Sacle Voltage
fin = 1431 * fs / num; % Input Signal Freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TI ADC Work Process %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% different channel sample
t         = zeros(n_ch, floor(num/n_ch)); % ts_sub = n_ch*ts
Vin_p_chs = zeros(n_ch, floor(num/n_ch));
Vin_n_chs = zeros(n_ch, floor(num/n_ch));
Vin_p_tot = zeros(1,num);
Vin_n_tot = zeros(1,num);

% ------ calculate total capacitor of stage1 SAR ADC ------- %
Cu1 = 20e-15; % ADC1 Unit Cap (F)
sigmaCu1 = sigmaC / sqrt(Cu1); % Unit Cap Mismatch
C_mismatch1 = sigmaCu1 * Cu1; % ADC1 Cap Mismatch (%)

C_arr_p1 = [2.^[(N1-1) :- 1:0], 1]; % Positive Plate Cap Array, CDAC_P   one more cap than stage2
C_arr_n1 = [2.^[(N1-1) :- 1:0], 1]; % Negative Plate Cap Array, CDAC_N

C_dev_p1 = C_mismatch1 * sqrt(C_arr_p1) .* randn(1,N1+1); % Positive Plate Cap Mismatch
C_dev_n1 = C_mismatch1 * sqrt(C_arr_n1) .* randn(1,N1+1); % Negative Plate Cap Mismatch

C_act_p1 = C_arr_p1 .* Cu1 + C_dev_p1; % Positive Plate Cap Array With Cap Mismatch
C_act_n1 = C_arr_n1 .* Cu1 + C_dev_n1; % Negative Plate Cap Array With Cap Mismatch

C_tot_p1 = sum(C_act_p1) + Cp_p1; % Positive Plate Total Cap With Parasitic Cap
C_tot_n1 = sum(C_act_n1) + Cp_n1; % Positive Plate Total Cap With Parasitic Cap
% ---------------------------------------------------------- %

% calibration（仅针对四通道）
delay = zeros(1, n_ch); % delay line
t0 = (0.05 * 1/fs) / 10; % the tuning step of the variable delay line
Vout_chs = zeros(n_ch + 1, floor(num/n_ch)); % 存储量化结果
Vout_diff = zeros(2, floor(num/n_ch)); % 量化结果差值通道
sign_delta = zeros(1, n_ch);
M = num/n_ch;
iteration = 40;
alpha0 = 1;
beta = 0;
delay_err_log = zeros(iteration, n_ch); % 每次迭代记录所有通道 delay 误差


for iter = 1 : iteration
    for i = 1 : n_ch
        % 采样时刻：t = ts * [i : n_ch : num]    e.g.[i:4:8]->[[1,5],[2,6],[3,7],[4,8]]
        t(i,:) = [ts*i : ts_sub : floor(num/n_ch)*ts_sub] + (Mis_TS(i)+delay(i)).*ones(1,floor(num/n_ch));
        t_sam = t(i,:);
        % 采样电压：
        Vin_p_chs(i,:) = Vref/2 + Mis_Gain(i).*(Vfs/2)*sin(2*pi*fin*t_sam) + Mis_OS(i) + sqrt(k*T/C_tot_p1)*randn(1,floor(num/n_ch));
        Vin_n_chs(i,:) = Vref/2 - Mis_Gain(i).*(Vfs/2)*sin(2*pi*fin*t_sam)             + sqrt(k*T/C_tot_n1)*randn(1,floor(num/n_ch));
    end
    % concat signals of different channels
    for i=1:n_ch
        for j = 1:floor(num/n_ch)
            Vin_p_tot(i + n_ch*(j - 1)) = Vin_p_chs(i,j);
        end
    end
    for i=1:n_ch
        for j = 1:floor(num/n_ch)
            Vin_n_tot(i + n_ch*(j - 1)) = Vin_n_chs(i,j);
        end
    end

    % Sub ADC
    Vresn = zeros(n_ch, floor(num/n_ch), N2+1);
    Vresp = zeros(n_ch, floor(num/n_ch), N2+1);
    D_chs = zeros(floor(num/n_ch), n_ch);
    for i = 1:n_ch
        Vin_n = Vin_n_chs(i,:);
        Vin_p = Vin_p_chs(i,:);

        [D,Vdacp,Vdacn] = pipelinedSAR_ADC(N1, ... %%%%%%%% ADC Params %%%%%%%%%
                            N2, ...
                            fs, ...
                            Vcm, ...
                            sigmaC, ...
                            Vref, ...            % Stage1 ADC Params
                            Cp_p1, ...
                            Cp_n1, ...
                            Vref, ...            % Stage2 ADC Params
                            Cp_p2, ...
                            Cp_n2, ...
                            k, ...               %%%%% Environment Params %%%%
                            T, ...
                            Av, ...              %%%%%%% ResAmp Params %%%%%%%
                            GBW, ...
                            Amp_offset, ...
                            Amp_noise, ...
                            Comp_offset1, ...    %%%%%% SAR Comp Params %%%%%%
                            Comp_offset2, ...
                            Comp_noise, ...
                            floor(num/n_ch),...  %%%%%%% Input Params %%%%%%%%
                            Vin_p, ...
                            Vin_n);
        % Store Results
        Vresn(i,:,:) = Vdacn;  
        Vresp(i,:,:) = Vdacp;
        D_chs(:,i)   = D;
        % 量化结果重建
        Vout_chs(i,:) = D .* Vref / 2^N;
    end
    Vout_chs(n_ch + 1, 1:end-1) = Vout_chs(1, 2:end);
    % 通过加法器和均值滤波器反馈到延时单元
    if iter < iteration / 2 % 第一阶段校正
        Vout_diff(1, 1:end-1) = abs( Vout_chs(1, 1:end-1) - Vout_chs(3, 1:end-1) );
        Vout_diff(2, 1:end-1) = abs( Vout_chs(3, 1:end-1) - Vout_chs(5, 1:end-1) );
        sign_delta(3) = -1*sign( sum( Vout_diff(1, 1:end-1)-Vout_diff(2, 1:end-1), "all" ) );
    elseif iter == iteration / 2
        sign_delta(3)=0; % 通道3校准已结束
    else % 第二阶段校正
        Vout_diff(1, 1:end-1) = abs( Vout_chs(1, 1:end-1) - Vout_chs(2, 1:end-1) );
        Vout_diff(2, 1:end-1) = abs( Vout_chs(2, 1:end-1) - Vout_chs(3, 1:end-1) );
        sign_delta(2) = -1*sign( sum( Vout_diff(1, 1:end-1)-Vout_diff(2, 1:end-1), "all" ) );
        Vout_diff(1, 1:end-1) = abs( Vout_chs(3, 1:end-1) - Vout_chs(4, 1:end-1) );
        Vout_diff(2, 1:end-1) = abs( Vout_chs(4, 1:end-1) - Vout_chs(5, 1:end-1) );
        sign_delta(4) = -1*sign( sum( Vout_diff(1, 1:end-1)-Vout_diff(2, 1:end-1), "all" ) );
    end

   
    
    % feedback, calibrate through delay line
    alpha = alpha0 / (1 + beta*iter);
    delay = delay + alpha*t0*sign_delta .* ones(1,n_ch);
    for i = 2:n_ch
    delay_err_log(iter, i) = delay(i) + Mis_TS(i) - Mis_TS(1);
    end

end

% 迭代曲线
figure;
hold on;
colors = lines(n_ch); % 不同颜色
for i = 1:n_ch
    plot(1:iteration, 1e12*delay_err_log(:, i), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', ['通道 ' num2str(i)]);
end
xlabel('迭代次数');
ylabel('delay(i) - Mis\_TS(i) / ps');
title('各通道校准残差变化（单位：ps）');
legend('show');
grid on;

% Mux
Dout = zeros(1,num);
for i = 1 : n_ch
    for j = 1 : floor(num / n_ch)
        Dout(i + n_ch * (j - 1)) = D_chs(j,i);
    end
end
Vout = Dout(1 : Nsample) .* Vref / 2^N;

% Plot Voltage Reconstructed From Output Codes
Show_Input = 1; % Show Input Signals? 0 : No, 1 : Yes
windowL_OUT = 0; % Left Boundary of Plot Window (us)
windowR_OUT = 2; % Right Boundary of Plot Window (us)
plot_Vout_Voltage(Vout, Vin_p_tot, Nsample, fs, windowL_OUT, windowR_OUT, Show_Input);

% Test and Plot Static & Dynamic Features
En_plot = 1; % Plot? ,0 : No, 1 : Yes
wid = 0; % Add Window? 0 : No, 1 : hamming, 2 : hann
[THD, SFDR, SNR, SNDR, ENOB] = Dynamic_test(Dout', fs, Nsample, En_plot, wid); % Dynamic Performance
[DNLmax, DNLmin, INLmax, INLmin] = Static_test(Dout', N); % Static Performance