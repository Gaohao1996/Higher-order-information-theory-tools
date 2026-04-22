function boxplot_trials_wins(M, varargin)
% 对“试次×窗口”矩阵 M 画箱线图

    p = inputParser;
    addParameter(p,'XLabels',[]);
    addParameter(p,'YLabel','Value',@ischar);
    addParameter(p,'Title','Boxplot by window',@ischar);
    parse(p, varargin{:});
    
    XLabels = p.Results.XLabels;
    YLabel  = p.Results.YLabel;
    TitleStr= p.Results.Title;

    [nTrials, nWins] = size(M);

    % 拉直
    Y = M(:);

    % 生成分组变量
    G = repelem(1:nWins, nTrials).';
    
    % 去掉 NaN
    valid = ~isnan(Y);
    Y = Y(valid);
    G = G(valid);

    % 👉 强制转成 categorical —— 这行能避免你遇到的所有报错
    G = categorical(G);

    % 画图
    figure('Color','w');
    boxplot(Y, G, 'Symbol','k+');
    grid on;

    ylabel(YLabel);
    title(TitleStr);

    if isempty(XLabels)
        XLabels = arrayfun(@(k) sprintf('W%d',k), 1:nWins, 'uni', 0);
    end
    set(gca,'XTick',1:nWins,'XTickLabel',XLabels,'FontSize',11);
end
