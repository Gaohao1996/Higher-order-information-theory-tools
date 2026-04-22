function [path_red, path_syn, g_red, g_syn, hist_red, hist_syn] = ...
    TE_syn_red(x, i, j, p, d, delayZ, phase)
% TE_syn_red
% ------------------------------------------------------------
% 对 i -> j 的传递熵进行冗余/协同路径分解
%
% 输入:
%   x     : [T x N] 原始数据矩阵
%   i     : driver 的原始编号，可以是标量或向量，例如 [1 2]
%   j     : target 的原始编号，可以是标量或向量，例如 [3 4]
%   p     : 模型阶数 / embedding order
%   d     : delay
%   phase : 1 表示相位数据, 0 表示普通连续数据
%
% 输出:
%   path_red  : 冗余路径中依次加入的条件变量原始编号
%   path_syn  : 协同路径中依次加入的条件变量原始编号
%   g_red     : 冗余路径上每一步 TE 值
%   g_syn     : 协同路径上每一步 TE 值
%   hist_red  : 冗余路径每一步的候选与选择历史
%   hist_syn  : 协同路径每一步的候选与选择历史
%
% 说明:
%   初始 TE:
%       g_red(1) = g_syn(1) = TE(i -> j)
%
%   之后每一步:
%       g_red(k+1) = min_c TE(i -> j | path_red, c)
%       g_syn(k+1) = max_c TE(i -> j | path_syn, c)
%
%   注意:
%   path_red/path_syn 中保存的是“条件路径”变量，而不是 driver
% ------------------------------------------------------------

    [~, nVar] = size(x);

    % 统一成行向量
    i = i(:)';
    j = j(:)';

    % -------- 输入合法性检查 --------
    if isempty(i)
        error('TE_syn_red: i cannot be empty.');
    end

    if isempty(j)
        error('TE_syn_red: j cannot be empty.');
    end

    if any(i < 1) || any(i > nVar)
        error('TE_syn_red: i contains out-of-range indices.');
    end

    if any(j < 1) || any(j > nVar)
        error('TE_syn_red: j contains out-of-range indices.');
    end

    if ~isempty(intersect(i, j))
        error('TE_syn_red: driver set i and target set j cannot overlap.');
    end

    % 固定的 driver / target 集合（原始编号）
    driver_set = i;
    target_set = j;

    % 候选条件变量（原始编号）
    candidate0 = setdiff(1:nVar, [driver_set target_set]);

    % ==========================================================
    % 初始 TE: TE(i -> j)
    % ==========================================================
    dt0 = [driver_set, target_set];
    x_sub0 = x(:, dt0);

    i_sub0 = 1:length(driver_set);
    j_sub0 = length(driver_set) + (1:length(target_set));

    TE0 = CTE_iteration(x_sub0, p, d, i_sub0, j_sub0, delayZ, phase);
    g0 = mean(TE0);

    g_red = g0;
    g_syn = g0;

    % 保存条件路径（原始编号）
    path_red = [];
    path_syn = [];

    % 各自独立的候选池
    ind_red = candidate0;
    ind_syn = candidate0;

    hist_red = struct([]);
    hist_syn = struct([]);

    % ==========================================================
    % 冗余路径：每次选择使 TE 最小的变量
    % ==========================================================
    step = 1;
    while ~isempty(ind_red)

        deltas = zeros(1, length(ind_red));

        for h = 1:length(ind_red)

            cand = ind_red(h);

            % 列顺序固定为:
            % [driver_set, target_set, 已选条件路径, 当前候选]
            dt = [driver_set, target_set, path_red, cand];
            x_sub = x(:, dt);

            % driver / target 的局部索引始终固定
            i_sub = 1:length(driver_set);
            j_sub = length(driver_set) + (1:length(target_set));

            TE_tmp = CTE_iteration(x_sub, p, d, i_sub, j_sub, delayZ, phase);
            deltas(h) = mean(TE_tmp);
        end

        [best_val, best_idx] = min(deltas);
        best_var = ind_red(best_idx);

        path_red(end+1) = best_var;
        g_red(end+1) = best_val;

        hist_red(step).step = step;
        hist_red(step).driver_set = driver_set;
        hist_red(step).target_set = target_set;
        hist_red(step).path_before = path_red(1:end-1);
        hist_red(step).candidates = ind_red;
        hist_red(step).delta_values = deltas;
        hist_red(step).chosen_var = best_var;
        hist_red(step).chosen_value = best_val;
        hist_red(step).path_after = path_red;

        ind_red(best_idx) = [];
        step = step + 1;
    end

    % ==========================================================
    % 协同路径：每次选择使 TE 最大的变量
    % ==========================================================
    step = 1;
    while ~isempty(ind_syn)

        deltas = zeros(1, length(ind_syn));

        for h = 1:length(ind_syn)

            cand = ind_syn(h);

            dt = [driver_set, target_set, path_syn, cand];
            x_sub = x(:, dt);

            i_sub = 1:length(driver_set);
            j_sub = length(driver_set) + (1:length(target_set));

            TE_tmp = CTE_iteration(x_sub, p, d, i_sub, j_sub, delayZ, phase);
            deltas(h) = mean(TE_tmp);
        end

        [best_val, best_idx] = max(deltas);
        best_var = ind_syn(best_idx);

        path_syn(end+1) = best_var;
        g_syn(end+1) = best_val;

        hist_syn(step).step = step;
        hist_syn(step).driver_set = driver_set;
        hist_syn(step).target_set = target_set;
        hist_syn(step).path_before = path_syn(1:end-1);
        hist_syn(step).candidates = ind_syn;
        hist_syn(step).delta_values = deltas;
        hist_syn(step).chosen_var = best_var;
        hist_syn(step).chosen_value = best_val;
        hist_syn(step).path_after = path_syn;

        ind_syn(best_idx) = [];
        step = step + 1;
    end
end









