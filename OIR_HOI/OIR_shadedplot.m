close all; clear; clc;

period = 'ict';    % 或 'ict'
load(sprintf('%s_humandata_10secwin.mat', period))

win_num = numel(OIR_results_human.result(1).win);
f = OIR_results_human.f(:);     % 列向量
index = [1 2 3 4 5 6 7];
N = numel(index);
M = numel(f);


%% frequecy band mask
% fmask = (f >= 2) & (f <= 50);

%% Plot
figure;
% T = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
colors = [31 119 180; 255 127 14; 214 39 40; 44 160 44]/255;
alpha_val = 0.18;
% 
for w = 1:win_num
    % —— 每个窗口重新收集该窗口的 7 个试次 —— 
    O12  = zeros(N, M);
    O1_2 = zeros(N, M);
    O2_1 = zeros(N, M);
    O1o2 = zeros(N, M);
    for k = 1:N
        i = index(k);
        O12(k,:)  = OIR_results_human.result(i).win(w).O12(:).';
        O1_2(k,:) = OIR_results_human.result(i).win(w).O1_2(:).';
        O2_1(k,:) = OIR_results_human.result(i).win(w).O2_1(:).';
        O1o2(k,:) = OIR_results_human.result(i).win(w).O1o2(:).';
    end
% 
%     nexttile; hold on;
%     % 如需全组件，解除下面三行注释
%     h1 = plot_iqr_shaded(f, O1_2,  colors(1,:), alpha_val);
%     h2 = plot_iqr_shaded(f, O1o2,  colors(2,:), alpha_val);
%     h3 = plot_iqr_shaded(f, O2_1,  colors(3,:), alpha_val);
%     h4 = plot_iqr_shaded(f, O12, colors(4,:), alpha_val);
%     
%     %% mask
% %     h1 = plot_iqr_shaded(f(fmask), O1_2(:,fmask),  colors(1,:), alpha_val);
% %     h2 = plot_iqr_shaded(f(fmask), O1o2(:,fmask),  colors(2,:), alpha_val);
% %     h3 = plot_iqr_shaded(f(fmask), O2_1(:,fmask),  colors(3,:), alpha_val);
% %     h4 = plot_iqr_shaded(f(fmask), O12(:,fmask), colors(4,:), alpha_val);
% 
% %     可视化核对：叠加 1 条单个试次看看是否与该窗口中位数一致
% %     plot(f, O12(1,:), 'k-', 'LineWidth', 0.7);  % 可选
%     title(sprintf('OIR density in %s period at window %d', period, w));
end
% 
% lgd = legend([h1 h2 h3 h4], ...
%     {'dir (sources→j)','coupling','reverse (j→sources)','total'}, ...
%     'Orientation','horizontal','NumColumns',4);
% 
% lgd.Layout.Tile = 'south';  % 放到整页底部

%% Whole time period
h1 = plot_iqr_shaded(f, O1_2,  colors(1,:), alpha_val);
h2 = plot_iqr_shaded(f, O1o2,  colors(2,:), alpha_val);
h3 = plot_iqr_shaded(f, O2_1,  colors(3,:), alpha_val);
h4 = plot_iqr_shaded(f, O12, colors(4,:), alpha_val);
xlabel('Frequency (Hz)');
ylabel('OIR spectral density (per Hz)');
title('OIR density in  ictal speriod in 10 sec');
lgd = legend([h1 h2 h3 h4], ...
    {'dir (sources→j)','coupling','reverse (j→sources)','total'}, ...
    'Orientation','horizontal','NumColumns',4);

% h1 = plot_iqr_shaded(f(fmask), O1_2(:,fmask),  colors(1,:), alpha_val);
% h2 = plot_iqr_shaded(f(fmask), O1o2(:,fmask),  colors(2,:), alpha_val);
% h3 = plot_iqr_shaded(f(fmask), O2_1(:,fmask),  colors(3,:), alpha_val);
% h4 = plot_iqr_shaded(f(fmask), O12(:,fmask), colors(4,:), alpha_val);
% title(sprintf('OIR density in %s period in 10 sec', period));
% lgd = legend([h1 h2 h3 h4], ...
%     {'dir (sources→j)','coupling','reverse (j→sources)','total'}, ...
%     'Orientation','horizontal','NumColumns',4);

