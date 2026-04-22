function Re = grouped_boxplot_resampled(mats, labels, K, varargin)
% 并排箱线图（重采样到统一 K 个相对时间段后）
% mats   : {M1, M2, ...}，每个 Mi 为 trials × windows 的矩阵（不同条件/窗口长度）
% labels : {'1s','2s',...} 条件标签（与 mats 对应）
% K      : 统一重采样到的相对时间段个数（比如 6 或 8）
%
% Name-Value 可选项：
%   'Interp'  : 'linear'|'nearest'|'pchip'（默认 'linear'）
%   'YLabel'  : y 轴标题（默认 'Value'）
%   'Title'   : 总标题（默认 'Grouped boxplots on relative time'）
%   'Test'    : 'none'|'ranksum'（仅两条件时可用；默认 'none'）
%   'FDR'     : true|false（对每个时间段的 p 做 BH-FDR；默认 true）
%   'Alpha'   : 显著性阈值（默认 0.05）
%
% 展示方式：每个相对时间段一个子图，子图内不同条件并排箱线图。

    p = inputParser;
    addParameter(p,'Interp','linear');
    addParameter(p,'YLabel','Value');
    addParameter(p,'Title','Grouped boxplots on relative time');
    addParameter(p,'Test','none');
    addParameter(p,'FDR',true);
    addParameter(p,'Alpha',0.05);
    parse(p,varargin{:});
    interpMethod = p.Results.Interp;
    ylab         = p.Results.YLabel;
    titleAll     = p.Results.Title;
    doTest       = p.Results.Test;
    doFDR        = p.Results.FDR;
    alpha        = p.Results.Alpha;

    nCond = numel(mats);
    assert(nCond == numel(labels), 'mats与labels数量不一致');

    %—— 重采样到 K ——%
    Re = cell(1, nCond); % 每个 trials × K
    for c = 1:nCond
        Re{c} = rebin_to_K(mats{c}, K, interpMethod);
    end

    %—— 统计检验（可选，仅两条件）——%
    pvals = nan(1,K);
    if strcmpi(doTest,'ranksum') && nCond == 2
        for k = 1:K
            v1 = Re{1}(:,k); v1 = v1(~isnan(v1));
            v2 = Re{2}(:,k); v2 = v2(~isnan(v2));
            if ~isempty(v1) && ~isempty(v2)
                pvals(k) = ranksum(v1, v2);
            end
        end
        if doFDR
            pvals_adj = bh_fdr(pvals);
        else
            pvals_adj = pvals;
        end
    else
        pvals_adj = [];
    end

    %—— 绘图：每个相对时间段一个子图，子图内并排箱线图 ——%
    figure('Color','w');
    for k = 1:K
        subplot(1,K,k);
        Yk = []; Gk = [];
        for c = 1:nCond
            Yk = [Yk; Re{c}(:,k)]; %#ok<AGROW>
            Gk = [Gk; repmat(c, size(Re{c},1), 1)]; %#ok<AGROW>
        end
        valid = ~isnan(Yk);
        Yk = Yk(valid);
        Gk = Gk(valid);

        boxplot(Yk, categorical(Gk), 'Symbol','k+'); grid on;
        set(gca,'XTickLabel',labels);
        if k == 1, ylabel(ylab); end
        titleTxt = sprintf('Relative window %d/%d', k, K);
        % 显著性标注（可选）
        if strcmpi(doTest,'ranksum') && nCond == 2 && ~isnan(pvals_adj(k))
            sig = stars(pvals_adj(k), alpha);
            title([titleTxt ' ' sig]);
        else
            title(titleTxt);
        end
    end
%     sgtitle(titleAll, 'FontSize', 10);
end

%——— 工具函数：行向量安全插值到统一 K 列 ———%
function R = rebin_to_K(M, K, method)
% M: trials × W（不同试次可含NaN）
% 结果：trials × K
    [nT, W] = size(M);
    if W == 0
        R = nan(nT, K); return;
    end
    x_src = linspace(0, 1, W);
    x_tar = linspace(0, 1, K);
    R = nan(nT, K);
    for t = 1:nT
        y = M(t,:);
        if all(isnan(y))
            continue;
        end
        % 先对NaN做一遍线性填补，避免插值失败（两端用邻值）
        y = fillmissing(y, 'linear', 'EndValues', 'nearest');
        R(t,:) = interp1(x_src, y, x_tar, method, 'extrap');
    end
end

%——— 工具函数：Benjamini–Hochberg FDR ———%
function p_adj = bh_fdr(p)
    p_adj = p;
    good = ~isnan(p);
    if ~any(good), return; end
    [ps, idx] = sort(p(good));
    m = numel(ps);
    q = ps .* m ./ (1:m);
    % 保序
    for i = m-1:-1:1
        q(i) = min(q(i), q(i+1));
    end
    tmp = nan(size(p));
    tmp(good) = 0;
    tmp(good) = q(invperm(idx));
    p_adj = tmp;
end

function y = invperm(idx)
% 返回排序索引的逆置换
    y = zeros(size(idx));
    y(idx) = 1:numel(idx);
end

%——— 工具函数：p值转星号 ———%
function s = stars(p, alpha)
    if nargin < 2, alpha = 0.05; end
    if isnan(p), s = ''; return; end
    if p < 0.001, s='***';
    elseif p < 0.01, s='**';
    elseif p < alpha, s='*';
    else, s='n.s.';
    end
end


