function plot_TE_path_significance(path_red, g_red, gc_red_surr, ...
                                   path_syn, g_syn, gc_syn_surr, ...
                                   source_variable, target_variable)

% ===================== Redundancy =====================
subplot(2,1,1); hold on;

Kred = length(path_red);
x_red = 0:Kred;                 % 第一个位置固定为 0
lab_red = [0, path_red];

%% surrogate 参考值（这里只画均值，你也可以改成中位数）
surr_red_mean = mean(gc_red_surr, 1);
% surr_red_mean = median(gc_red_surr, 1);

% 冗余显著性：真实值显著小于 surrogate
sig_red = false(1, Kred+1);
sig_red(1) = true;   % 起点总是高亮显示

p_red = nan(1, Kred);
for dr = 1:Kred
    p_red(dr) = mean(gc_red_surr(:,dr) <= g_red(dr+1));
    sig_red(dr+1) = p_red(dr) < 0.05;
end

% 先画 surrogate（空心圈，紫色）
scatter(x_red(2:end), surr_red_mean, 80, 'o', ...
    'MarkerEdgeColor', [0 1 0], 'LineWidth', 1.5);

% 再画真实点
for k = 1:(Kred+1)
    if sig_red(k)
        scatter(x_red(k), g_red(k), 60, 'r', 'filled');
    else
        scatter(x_red(k), g_red(k), 60, 'k', 'filled');
    end
end

xticks(x_red);
xticklabels(lab_red);
ylabel('TM (redundancy)');
xlabel(['driver ', num2str(source_variable), ', target ', num2str(target_variable)]);
title('Redundancy path');
box on;


% ===================== Synergy =====================
subplot(2,1,2); hold on;

Ksyn = length(path_syn);
x_syn = 0:Ksyn;
lab_syn = [0, path_syn];

% surr_syn_mean = mean(gc_syn_surr, 1);
surr_syn_mean = median(gc_syn_surr, 1);

% 协同显著性：真实值显著大于 surrogate
sig_syn = false(1, Ksyn+1);
sig_syn(1) = true;

p_syn = nan(1, Ksyn);
for dr = 1:Ksyn
    p_syn(dr) = mean(gc_syn_surr(:,dr) >= g_syn(dr+1));
    sig_syn(dr+1) = p_syn(dr) < 0.05;
end

% surrogate（空心圈，紫色）
scatter(x_syn(2:end), surr_syn_mean, 80, 'o', ...
    'MarkerEdgeColor', [0 1 0], 'LineWidth', 1.5);

% 真实点
for k = 1:(Ksyn+1)
    if sig_syn(k)
        scatter(x_syn(k), g_syn(k), 60, 'b', 'filled');
    else
        scatter(x_syn(k), g_syn(k), 60, 'k', 'filled');
    end
end

xticks(x_syn);
xticklabels(lab_syn);
ylabel('TM (synergy)');
xlabel(['driver ', num2str(source_variable), ', target ', num2str(target_variable)]);
title('Synergy path');
box on;

end
