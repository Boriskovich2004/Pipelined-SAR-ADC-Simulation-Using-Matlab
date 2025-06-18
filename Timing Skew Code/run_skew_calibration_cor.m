% 基于自相关的时钟校正
function [SNDR_out, ENOB_out] = run_skew_calibration_cor(iteration, seed)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% ADC Parameter %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N1 = 6;          % stage1 SAR
    N2 = 8;          % stage2 SAR
    N = N1 + N2 - 1; % resolution (one bit redundancy)
    fs = 100e6; % Sample Rate(MHz)
    ts = 1/(fs);
    Vref = 1.8; % Reference Voltage
    Vcm = Vref/2; % Common Mode Voltage
    LSB = Vref/(2^N);
    n_ch = 4; % number of channel

    %%%%%%%%%% Mismatch & noise %%%%%%%%%%%
    % Inter-Channel Mismatch
    Mis_OS = 000 * LSB * randn(1,n_ch);       % Inter-channel Offset
    Mis_Gain = 1 + 000 * LSB * randn(1,n_ch); % Inter-channel Gain Mismatch
    rng(seed);
    Mis_TS = 0.005 * 1/fs *  randn(1, n_ch);    % Inter-channel Clock Skew

    % SAR Capacitor Mismatch
    sigmaC = 0.00;    % Global Capacitor Mismatch
    Cp_p1  = 000e-15; % Stage1 Positive Plate Parasitic Capacitor (F)
    Cp_n1  = 000e-15; % Stage1 Negative Plate Parasitic Capacitor (F)
    Cp_p2  = 000e-15; % Stage2 Positive Plate Parasitic Capacitor (F)
    Cp_n2  = 000e-15; % Stage2 Negative Plate Parasitic Capacitor (F)

    % ResAmp Mismatch & noise
    Av         = 120;    % Open-loop Gain (dB)
    GBW        = 10e9;   % Gain-Bandwidth Product (GHz)
    Amp_offset = 000e-3; % Amplifier Offset (V)
    Amp_noise  = 000e-3; % Amplifier RMS Voltage Noise (V)

    % SAR Comparator Mismatch & noise
    Comp_offset1 = 000e-3;% Stage1 SAR Comp Offset (V)
    Comp_offset2 = 000e-3;% Stage2 SAR Comp Offset (V)
    Comp_noise   = 000e-3;% SAR Comp RMS Voltage Noise (V)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Environment Parameter %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = 1.38e-23;       % Bolzman Constant
    T = 000;            % Temperature (K)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Sub ADC Parameter %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fs_sub = fs / n_ch;    % Sub ADC Sample Rate
    ts_sub = 1 / (fs_sub); % Sub ADC Sample Cycle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% input signal %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num = 2^15;            % Sample Points
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
    Vin_p_tot = zeros(num);
    Vin_n_tot = zeros(num);

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
    delay = zeros(1, n_ch);
    t0 = (0.01 * 1/fs) / 100; % the tuning step of the variable delay line
    Vout_chs = zeros(n_ch, floor(num/n_ch));
    R_ti = zeros(1, n_ch);
    R_ts = 0;
    R_delta = zeros(1, n_ch);
    sign_delta = zeros(1, n_ch);
    M = num/n_ch;
    iteraion = 200;
    alpha0 = 1;
    beta = 0;
    delay_err_log = zeros(iteraion, n_ch); % 每次迭代记录所有通道 delay 误差


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
        % 计算相邻样本的自相关函数
        for i =1:n_ch
            if i<n_ch
                R_ti(i) = sum( Vout_chs(i,1:M-1) .* Vout_chs(i+1,1:M-1) ) / (M-1);
            else
                R_ti(i) = sum( Vout_chs(n_ch,1:M-1) .* Vout_chs(1,2:M) ) / (M-1);
            end
        end
        % calculate calibration orientation 
        R_ts = sum(R_ti, "all")/n_ch;
        for i=1:n_ch
        R_delta(i) = sum(R_ti(1:i-1)) - (i-1)*R_ts;
        end
        sign_delta = sign(R_delta);
        % feedback, calibrate through delay line
        alpha = alpha0 / (1 + beta*iter);
        delay = delay + alpha*t0*sign_delta .* ones(1,n_ch);
        for i = 1:n_ch
        delay_err_log(iter, i) = delay(i) + Mis_TS(i) - Mis_TS(1);
        end

    end

    % Mux
    Dout = zeros(1,num);
    for i = 1 : n_ch
        for j = 1 : floor(num / n_ch)
            Dout(i + n_ch * (j - 1)) = D_chs(j,i);
        end
    end
    Vout = Dout(1 : Nsample) .* Vref / 2^N;


    % Test and Plot Static & Dynamic Features
    En_plot = 1; % Plot? ,0 : No, 1 : Yes
    wid = 0; % Add Window? 0 : No, 1 : hamming, 2 : hann
    [~, ~, ~, SNDR_out, ENOB_out] = Dynamic_test(Dout', fs, Nsample, En_plot, wid); % Dynamic Performance
    [DNLmax, DNLmin, INLmax, INLmin] = Static_test(Dout', N); % Static Performance


end