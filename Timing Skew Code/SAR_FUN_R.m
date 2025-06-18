% clc;
% clear;
% close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% ADC基本参数 %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N = 10; % ADC分辨率(Bit)
% fs = 100e6; % ADC采样率(Hz)
% Vcm = 0.9; %共模电压(V)
% sigmaC = 0.00; % 电容适配
%
% %%%%%%%%%%%%% 定义第一级ADC参数 %%%%%%%%%%%%%
% Cu = 20e-15; % ADC1 单位电容(F)
% sigmaCu = sigmaC/sqrt(20);
% C_mismatch = sigmaCu * Cu; % ADC1 电容失配(%)
% C_arr_p = [2.^[(N-2) :- 1:0],1]; % 电容阵列正极板, CDAC_P
% C_arr_n = [2.^[(N-2) :- 1:0],1]; % 电容阵列负极板, CDAC_N
% C_dev_p = C_mismatch*sqrt(C_arr_p) .* randn(1,N); % 正极板电容失配 比第二级多一个电容
% C_dev_n = C_mismatch*sqrt(C_arr_n) .* randn(1,N); % 负极板电容失配
% C_act_p = C_arr_p .* Cu+C_dev_p; % 电容阵列正极板,考虑电容失配
% C_act_n = C_arr_n .* Cu+C_dev_n; % 电容阵列负极板,考虑电容失配
% Cp_p1 = 0e-15; % 正极板寄生电容(F)
% Cp_n1 = 0e-15; % 负极板电容失配(F)
% C_tot_p = sum(C_act_p)+Cp_p1; % 总电容大小,考虑寄生电容
% C_tot_n = sum(C_act_n)+Cp_n1; % 总电容大小,考虑寄生电容
% Vref=1.8; % ADC1参考电压(V)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% 定义环境参数 %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% k = 1.38e-23; % 玻尔兹曼常数
% Jitter = 000e-15; % 时钟抖动(s)
% T = 300; % 温度(K)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%% 定义输入信号参数 %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% num = 2^15; % 采样点数
% Vfs = 1.8; % 输入信号摆幅(V)
% % fin = 1e6; % 输入信号频率(Hz)
% fin = 1411/num*fs; % 输入信号频率(Hz)
% ts = 1/(fs); % 输入信号采样周期(s)
% t = [0:ts:(num-1)*ts]'+Jitter .* randn(num,1); % 采样序列
% Vin_p = Vcm+(Vfs/2)*sin(2*pi*fin*t); % 正输入端采样值
% Vin_n = Vcm-(Vfs/2)*sin(2*pi*fin*t); % 负输入端采样值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% ADC的工作过程 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D,Vdacp,Vdacn] = SAR_FUN_R(N, num, Vref, Vin_p, Vin_n, k, T, C_tot_p, C_tot_n, C_act_p, C_act_n, Comp_offset, Comp_noise)
Vresp = zeros(1,N+1);
Vresn = zeros(1,N+1);
Vdacp = zeros(num,N+1);
Vdacn = zeros(num,N+1);
B = zeros(1,N);
D = zeros(num,N);
for j=1:num %使用for循环来实现周期工作过程
    Vip=Vin_p(j);%正输入端电压,单次
    Vin=Vin_n(j);%负输入端电压,单次
    % if(j==2048)
    %     test=0;
    % end
    Vresp(1)=Vip;
    Vresn(1)=Vin;
    for i=1:N % N+1次比较，比最后一级多一次
        if(Vip-Vin-Comp_offset-Comp_noise*randn() <= 0)
            B(i)=0;
            Vip=Vip+Vref/2*C_act_p(i)/C_tot_p; % Vcm-BasedlFf
            Vin=Vin-Vref/2*C_act_n(i)/C_tot_n;
            % Vip=Vip+Vref*C_act_p(i)/C_tot_p; % MonotonicHFf
            % % Vin=Vin-Vref*C_act_n(i)/C_tot_n;
            Vresp(i+1)=Vip;
            Vresn(i+1)=Vin;
        else
            B(i)=1;
            Vip=Vip-Vref/2*C_act_p(i)/C_tot_p; % Vcm-BasedBfFf
            Vin=Vin+Vref/2*C_act_n(i)/C_tot_n;
            % Vip=Vip-Vref*C_act_p(i)/C_tot_p; % MonotonicH/f
            % % Vin=Vin+Vref/2*C_act_n(i)/C_tot_n;
            Vresp(i+1)=Vip;
            Vresn(i+1)=Vin;
        end
    end
    % if (Vip-Vin-Comp_offset-Comp_noise*randn() <= 0)
    %     B(N)=0;
    %     Vresp(N+1)=Vip;
    %     Vresn(N+1)=Vin;
    % else
    %     B(N)=1;
    %     Vresp(N+1)=Vip;
    %     Vresn(N+1)=Vin;
    %
    % end

    D(j,:)=B;% ADC的数字码输出
    Vdacp(j,:)=Vresp;%ADC的正端CBAC上的余差电压
    Vdacn(j,:)=Vresn;%ADC的负端CDAC上的余差电压
end

% Dout=D*2.^[(N-1) :- 1:0]';%数字码经过理想DAC复原之后的输出
% Vout=Dout*Vref/2^N;%归一化之后的输出

end

% %%%%%%%%% 测试 %%%%%%%%%%%
% Nsample = num;% FFT点数
%
% % 绘制CDAC余差电压波形
% Show_Input = 0; % 是否展示输入信号,0=No,1=Yes
% windowL_DAC = 0; % 视窗左边界(us)
% windowR_DAC = 2; % 视窗右边界(us)
% [vector_Vresp, vector_Vresn] = plot_CDAC_Voltage(Vdacp, Vdacn, Vin_p, Vin_n, N, Nsample, fs, windowL_DAC, windowR_DAC, Show_Input);
%
% % 绘制输出复原后的电压波形
% Show_Input = 1; % 是否展示输入信号,0=No,1=Yes
% windowL_OUT = 0; % 视窗左边界(us)
% windowR_OUT = 2; % 视窗右边界(us)
% plot_Vout_Voltage(Vout, Vin_p, Nsample, fs, windowL_OUT, windowR_OUT, Show_Input);
%
% % 测试并绘制动态和静态特性
% En_plot = 1; % 是否绘图,0=No,1=Yes
% wid = 0; % 是否加窗,0=No,1=hamming, 2=hann
% [THD, SFDR, SNR, SNDR, ENOB] = Dynamic_test(Dout', fs, Nsample, En_plot, wid); % 动态性能测试
% [DNLmax, DNLmin, INLmax, INLmin] = Static_test(Dout', N); % 静态性能测试
