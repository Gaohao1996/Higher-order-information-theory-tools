%% === 输入：U, R, S 为 7x4 矩阵（行=试次，列=频段） ===
% 例：已在工作区有 U, R, S
U = load('U_total_pre.mat').U_information_total;
R = load('R_total_pre.mat').R_information_total;
S = load('S_total_pre.mat').S_information_total;

% ---------- 1) 逐(试次×频段)归一化 ----------
den = U + R + S;                    % 7x4
den(den == 0) = NaN;                % 防零除；若真的为0，则该格会被忽略
U_n = U ./ den;
R_n = R ./ den;
S_n = S ./ den;


figure; hold on;

% === 1) 画 boxplot：关闭所有散点（Outliers）===
boxplot(U_n, 'Positions', posU, 'Widths', w, ...
    'Colors',[0 0 0], 'Symbol','');   % Unique → 黑

boxplot(R_n, 'Positions', posR, 'Widths', w, ...
    'Colors',[1 0 0], 'Symbol','');   % Redundant → 红

boxplot(S_n, 'Positions', posS, 'Widths', w, ...
    'Colors',[0 0 1], 'Symbol','');   % Synergistic → 蓝


% === 2) 统一线宽（可选，但好看）===
set(findobj(gca,'Tag','Median'),'LineWidth',1.8);
set(findobj(gca,'Tag','Box'),'LineWidth',1.2);
set(findobj(gca,'Tag','Whisker'),'LineWidth',1.2);


% === 3) 坐标轴与标签 ===
xlim([0.5, numel(freq_names)+0.5]);
set(gca,'XTick',pos_base,'XTickLabel',freq_names,'FontSize',12);
ylabel('Normalized proportion (U / R / S)','FontSize',12);
title('Transfer entropy decomposition in preictal period (normalized)','FontSize',13);


% === 4) 彻底不用自动 legend，只手动画 3 条 ===
legend off

hU = plot(nan,nan,'-','LineWidth',4,'Color',[0 0 0]);
hR = plot(nan,nan,'-','LineWidth',4,'Color',[1 0 0]);
hS = plot(nan,nan,'-','LineWidth',4,'Color',[0 0 1]);

legend([hU hR hS], ...
    {'U (Unique)','R (Redundant)','S (Synergistic)'}, ...
    'Location','northeast');

box on; grid on;
