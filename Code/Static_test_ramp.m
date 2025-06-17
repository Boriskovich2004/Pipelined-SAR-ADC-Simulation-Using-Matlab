function [DNLmax, DNLmin, INLmax, INLmin] = Static_test_ramp(code, Nbit)
    % 使用理想斜坡输入进行静态测试计算 DNL / INL 并保存图像
    
    % 限制输出码范围
    code = round(code);
    code(code < 0) = 0;
    code(code > 2 ^ Nbit - 1) = 2 ^ Nbit - 1;

    % 获取待测 ADC 的直方图
    % [counts, ~] = histcounts(code, 0:2 ^ Nbit);
    [counts, ~] = histcounts(code, -0.5:1:(2^Nbit - 0.5));
    M = sum(counts); % 总采样点数

    % === 使用斜坡输入生成 ideal_counts === %
    ideal_counts = ones(1, 2 ^ Nbit) * (M / 2 ^ Nbit);

    % === DNL 和 INL 计算 === %
    DNL = (counts - ideal_counts) ./ ideal_counts;
    INL = cumsum(DNL);
    % 
    % % === 强制 INL 首尾为 0 的线性趋势移除 === %
    % N = length(INL_raw);
    % linear_trend = linspace(INL_raw(1), INL_raw(end), N);  % 从 INL(1) 到 INL(end) 的直线
    % INL = INL_raw - linear_trend;  

    fprintf('INL at code 0     = %.2f LSB\n', INL(1));
    fprintf('INL at code end   = %.2f LSB\n', INL(end));

    % 最大/最小值输出
    DNLmax = max(DNL);
    DNLmin = min(DNL);
    INLmax = max(INL);
    INLmin = min(INL);

    % === 可视化 === %
    fig = figure('Visible', 'on'); % 可设为 'off' 以隐藏弹窗

    subplot(2, 1, 1);
    bar(DNL);
    xlabel('Code');
    ylabel('DNL (LSB)');
    title('DNL (Ramp Test)');
    grid on;
    txt1 = {
            sprintf('DNLmax   = %.2f LSB', DNLmax)
            sprintf('DNLmin   = %.2f LSB', DNLmin)
            };
    text(length(DNL) * 0.8, max(DNL) * 0.7, txt1, 'FontSize', 5, 'BackgroundColor', 'w', 'EdgeColor', 'k');

    subplot(2, 1, 2);
    bar(INL);
    xlabel('Code');
    ylabel('INL (LSB)');
    title('INL (Cumulative, Ramp Test)');
    grid on;
    txt2 = {
            sprintf('INLmax   = %.2f LSB', INLmax)
            sprintf('INLmin   = %.2f LSB', INLmin)
            };
    text(length(INL) * 0.8, max(INL) * 0.9, txt2, 'FontSize', 5, 'BackgroundColor', 'w', 'EdgeColor', 'k');

    % % === 保存图像 === %
    % outputFolder = 'Pictures';
    % 
    % if ~exist(outputFolder, 'dir')
    %     mkdir(outputFolder);
    % end
    % 
    % timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    % filename = fullfile(outputFolder, ['StaticTest_' timestamp '.png']);
    % saveas(fig, filename);
    % 
    % disp(['✅ 图像已保存到: ' filename]);

end
