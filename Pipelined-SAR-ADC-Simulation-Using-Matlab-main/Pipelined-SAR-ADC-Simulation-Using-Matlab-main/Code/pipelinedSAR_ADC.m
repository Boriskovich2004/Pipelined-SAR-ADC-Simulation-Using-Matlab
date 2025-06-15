% clc;
% clear;
% close all;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%% ADC Parameter %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N1 = 6;        % stage1 SAR
% N2 = 8;        % stage2 SAR
% fs = 100e6;    % Sample Rate(MHz)
% Vcm = 0.9;     % Common Mode Voltage
% sigmaC = 0.00; % Global Capacitor Mismatch
%
% %%%%%%%%% Stage1 ADC Parameter %%%%%%%%%%
% Vref1=1.8; % ADC1 Reference Voltage (V)
% Cp_p1 = 0e-15; % Stage1 Positive Plate Parasitic Capacitor(F)
% Cp_n1 = 0e-15; % Stage1 Negative Plate Parasitic Capacitor(F)
%
% %%%%%%%%% Stage2 ADC Parameter %%%%%%%%%%
% Vref2 = 1.8; % ADC2 Reference Voltage (V)
% Cp_p2 = 0e-15; % Positive Plate Parasitic Capacitor (F)
% Cp_n2 = 0e-15; % Negative Plate Parasitic Capacitor (F)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Environment Parameter %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k = 1.38e-23;     % Bolzman Constant
% T = 000;          % Temperature (K)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% Amplifier Parameter %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Av         = 120;    % Open-loop Gain (dB)
% GBW        = 10e9;   % Gain-Bandwidth Product (GHz)
% Amp_offset = 000e-3; % Amplifier Offset (V)
% Amp_noise  = 000e-3; % Amplifier RMS Voltage Noise (V)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% Comparater Parameter %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comp_offset1 = 000e-3;% Stage1 SAR Comp Offset (V)
% Comp_offset2 = 000e-3;% Stage2 SAR Comp Offset (V)
% Comp_noise   = 000e-3;% SAR Comp RMS Voltage Noise (V)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% input signal %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num = 2^15;           % Sample Points
% Vfs = 1.8;            % Input Signal Full Sacle Voltage (V)
% % fin = 1e6;          % Input Signal Freq (Hz)
% fin = 1511/num*fs;    % Input Signal Freq (Hz)
% ts = 1/(fs);          % Input Signal Cycle (s)
% Jitter = 000e-15;     % Clock Jitter (s)
% t = [0:ts:(num-1)*ts]'+Jitter .* randn(num,1);    % Sample Sequence
% Vin_p = Vcm+(Vfs/2)*sin(2*pi*fin*t);              % Positive Input
% Vin_n = Vcm-(Vfs/2)*sin(2*pi*fin*t);              % Negative Input
% % Vin_p = Vcm+(Vfs/2)*(2*t/((num-1)*ts)-1);       % Positive Input
% % Vin_n = Vcm-(Vfs/2)*(2*t/((num-1)*ts)-1);       % Negative Input
%
% [D,Vdacp,Vdacn] = SAR(N1, ...            %%%%%%%% ADC Params %%%%%%%%%
%                       N2, ...
%                       fs, ...
%                       Vcm, ...
%                       sigmaC, ...
%                       Vref1, ...         % Stage1 ADC Params
%                       Cp_p1, ...
%                       Cp_n1, ...
%                       Vref2, ...         % Stage2 ADC Params
%                       Cp_p2, ...
%                       Cp_n2, ...
%                       k, ...             %%%%% Environment Params %%%%
%                       T, ...
%                       Av, ...            %%%%%%% ResAmp Params %%%%%%%
%                       GBW, ...
%                       Amp_offset, ...
%                       Amp_noise, ...
%                       Comp_offset1, ...  %%%%%% SAR Comp Params %%%%%%
%                       Comp_offset2, ...
%                       Comp_noise, ...
%                       num, ...           %%%%%%% Input Params %%%%%%%%
%                       Vin_p, ...
%                       Vin_n);
%
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% ADC Test & Plot %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N  = N1 + N2 - 1;       % Total Resolution (1 bit Redundancy)
% Vout = D * Vref1 / 2^N; % Normalized Output
% Nsample = num;          % FFT Points
%
% % Plot CDAC Residue Voltage
% Show_Input = 0; % Show Input? 0 : No, 1 : Yes
% windowL_DAC = 0; % Left Boundary of Plot Window (us)
% windowR_DAC = 2; % Right Boundary of Plot Window (us)
% [vector_Vresp, vector_Vresn] = plot_CDAC_Voltage(Vdacp, Vdacn, Vin_p, Vin_n, N2, Nsample, fs, windowL_DAC, windowR_DAC, Show_Input);
%
% % Plot Voltage Reconstructed From Output Codes
% Show_Input = 1; % Show Input Signals? 0 : No, 1 : Yes
% windowL_OUT = 0; % Left Boundary of Plot Window (us)
% windowR_OUT = 2; % Right Boundary of Plot Window (us)
% plot_Vout_Voltage(Vout, Vin_p, Nsample, fs, windowL_OUT, windowR_OUT, Show_Input);
%
% % Test and Plot Static & Dynamic Features
% En_plot = 1; % Plot? ,0 : No, 1 : Yes
% wid = 0; % Add Window? 0 : No, 1 : hamming, 2 : hann
% [THD, SFDR, SNR, SNDR, ENOB] = Dynamic_test(D', fs, Nsample, En_plot, wid); % Dynamic Performance
% [DNLmax, DNLmin, INLmax, INLmin] = Static_test(D', N); % Static Performance
% function [D,Vdacp,Vdacn] = SAR(N1, ...            %%%%%%%% ADC Params %%%%%%%%%
%                       N2, ...
%                       fs, ...
%                       Vcm, ...
%                       sigmaC, ...
%                       Vref1, ...         % Stage1 ADC Params
%                       Cp_p1, ...
%                       Cp_n1, ...
%                       Vref2, ...         % Stage2 ADC Params
%                       Cp_p2, ...
%                       Cp_n2, ...
%                       k, ...             %%%%% Environment Params %%%%
%                       T, ...
%                       Av, ...            %%%%%%% ResAmp Params %%%%%%%
%                       GBW, ...
%                       Amp_offset, ...
%                       Amp_noise, ...
%                       Comp_offset1, ...  %%%%%% SAR Comp Params %%%%%%
%                       Comp_offset2, ...
%                       Comp_noise, ...
%                       num, ...           %%%%%%% Input Params %%%%%%%%
%                       Vin_p, ...
%                       Vin_n)

function [D,Vdacp,Vdacn] = pipelinedSAR_ADC(N1, ...            %%%%%%%% ADC Params %%%%%%%%%
    N2, ...
    fs, ...
    Vcm, ...
    sigmaC, ...
    Vref1, ...         % Stage1 ADC Params
    Cp_p1, ...
    Cp_n1, ...
    Vref2, ...         % Stage2 ADC Params
    Cp_p2, ...
    Cp_n2, ...
    k, ...             %%%%% Environment Params %%%%
    T, ...
    Av, ...            %%%%%%% ResAmp Params %%%%%%%
    GBW, ...
    Amp_offset, ...
    Amp_noise, ...
    Comp_offset1, ...  %%%%%% SAR Comp Params %%%%%%
    Comp_offset2, ...
    Comp_noise, ...
    num, ...           %%%%%%% Input Params %%%%%%%%
    Vin_p, ...
    Vin_n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ADC Parameter %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = N1 + N2 - 1; % ADC Resolution (bit)

%%%%%%%%% Stage1 ADC Parameter %%%%%%%%%%
Cu1 = 20e-15; % ADC1 Unit Cap (F)
sigmaCu1 = sigmaC / sqrt(Cu1); % ADC1 Unit Cap Mismatch
C_mismatch1 = sigmaCu1 * Cu1; % ADC1 Cap Mismatch (%)

C_arr_p1 = [2.^[(N1-1) :- 1:0], 1]; % Positive Plate Cap Array, CDAC_P   one more cap than stage2
C_arr_n1 = [2.^[(N1-1) :- 1:0], 1]; % Negative Plate Cap Array, CDAC_N

C_dev_p1 = C_mismatch1 * sqrt(C_arr_p1) .* randn(1,N1+1); % Positive Plate Cap Mismatch
C_dev_n1 = C_mismatch1 * sqrt(C_arr_n1) .* randn(1,N1+1); % Negative Plate Cap Mismatch

C_act_p1 = C_arr_p1 .* Cu1 + C_dev_p1; % Positive Plate Cap Array With Cap Mismatch
C_act_n1 = C_arr_n1 .* Cu1 + C_dev_n1; % Negative Plate Cap Array With Cap Mismatch

C_tot_p1 = sum(C_act_p1) + Cp_p1; % Positive Plate Total Cap With Parasitic Cap
C_tot_n1 = sum(C_act_n1) + Cp_n1; % Positive Plate Total Cap With Parasitic Cap

%%%%%%%%% Stage2 ADC Parameter %%%%%%%%%%
Cu2 = 1e-15; % ADC1 Unit Cap (F)
sigmaCu2 = sigmaC; % ADC2 Unit Cap Mismatch
C_mismatch2 = sigmaCu2*Cu2; % ADC2 Cap Mismatch (%)

C_arr_p2 = [2.^[(N2-2) :- 1:0],1]; % Positive Plate Cap Array, CDAC_P   one less cap than stage1
C_arr_n2 = [2.^[(N2-2) :- 1:0],1]; % Negative Plate Cap Array, CDAC_N

C_dev_p2 = C_mismatch2*sqrt(C_arr_p2) .* randn(1,N2); % Positive Plate Cap Mismatch
C_dev_n2 = C_mismatch2*sqrt(C_arr_n2) .* randn(1,N2); % Negative Plate Cap Mismatch

C_act_p2 = C_arr_p2 .* Cu2+C_dev_p2; % Positive Plate Cap Array With Cap Mismatch
C_act_n2 = C_arr_n2 .* Cu2+C_dev_n2; % Negative Plate Cap Array With Cap Mismatch

C_tot_p2 = sum(C_act_p2)+Cp_p2; % Positive Plate Total Cap With Parasitic Cap
C_tot_n2 = sum(C_act_n2)+Cp_n2; % Positive Plate Total Cap With Parasitic Cap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Pipelined SAR ADC Work Process %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stage1 SAR ADC
[D1,Vdacp1,Vdacn1] = SAR_FUN_R(N1, ...
    num, ...
    Vref1, ...
    Vin_p, ...
    Vin_n, ...
    k, ...
    T, ...
    C_tot_p1, ...
    C_tot_n1, ...
    C_act_p1, ...
    C_act_n1, ...
    Comp_offset1, ...
    Comp_noise);
% Store Residue
Vresp1 = Vdacp1(:,N1+1)';
Vresn1 = Vdacn1(:,N1+1)';

% ResAmp
Cs1 = 2 * Cu1 * Vref1 / Vref2; % Sample Cap
Ch1 = 2^N1 * Cu1 + Cp_n1; % Hold Cap With Parasitic Cap
% FeedBack Gain = Cs1 / Ch1, Close-loop Gain ~ Ch1 / Cs1 ~ 2^N1
Vinp1 = Vresp1;
Vinn1 = Vresn1;
[Voutp1,Voutn1] = RA_FUN(Vinp1, ...
    Vinn1, ...
    Av, ...
    GBW, ...
    Cs1, ...
    Ch1, ...
    Vcm, ...
    Amp_offset, ...
    Amp_noise, ...
    fs);

% Stage2 SAR ADC
[D2,Vdacp,Vdacn] = SAR_FUN(N2, ...
    num, ...
    Vref2, ...
    Voutp1, ...
    Voutn1, ...
    k, ...
    T, ...
    C_tot_p2, ...
    C_tot_n2, ...
    C_act_p2, ...
    C_act_n2, ...
    Comp_offset2, ...
    Comp_noise);

% Data Processing
D_d1 = D1*2.^[(N - 1) :- 1:(N - N1)]'; % Stage1 Output in Decimal
D_d2 = D2*2.^[(N - N1) :- 1:(N - N1 - N2 + 1)]'; % Stage2 Output in Decimal
D = D_d1 + D_d2 - 2.^(N - N1 - 1); % Total Result with Redundancy

% figure(1)
% plot(t, DOUT1*Vref1/2^N, 'b');
% hold on;
% plot(t, DOUT2*Vref1/2^N, 'r');
% plot(t, Vin_p, 'k');
% plot(t,Voutp1, 'c');
end