function [ax, h] = boxplot_diy(M, varargin)
% BOX PLOT DIY — trials × windows 的箱线图（支持缩放与自动刻度）
%
% 用法示例：
%   [ax,h] = boxplot_diy(M, 'Parent', subplot(2,2,1), ...
%       'XLabels', {'W1','W2','W3','W4','W5'}, ...
%       'YLabel', 'O-information (per Hz)', ...
%       'Title',  'Theta band', ...
%       'YScaleFactor', 1e3, ...        % 仅用于显示放大（坐标轴会标注 ×10^{-3}）
%       'AutoTight', true, ...
%       'Colors', [0.20 0.60 0.20], ... % 或 nWins×3
%       'BoxFaceAlpha', 0.25, ...
%       'Symbol','k+');
%
% 输入：
%   M : [nTrials × nWins]，每列对应一个窗口（W1..Wn）
%
% 关键参数（Name-Value）：
%   'Parent'        : 目标 axes 句柄，默认 gca
%   'XLabels'       : 每个窗口的标签，默认 {'W1','W2',...}
%   'YLabel'        : y 轴标题（字符串），默认 'Value'
%   'Title'         : 图标题（字符串），默认 'Boxplot by window'
%   'Colors'        : 1×3 或 nWins×3 的 RGB（0~1），用于边框与填充
%   'BoxFaceAlpha'  : 填充透明度 0~1（>0 才会填充）
%   'Symbol'        : 离群点标记（传给 boxplot）
%   'YScaleFactor'  : 仅用于显示的纵轴缩放（默认 1，不改数据）
%   'YLim'          : [ymin ymax] 手动设定 y 轴范围（作用于显示后的尺度）
%   'AutoTight'     : true/false，自动根据数据贴边设置 y 轴范围（默认 true）
%
% 返回：
%   ax : 使用的坐标轴句柄
%   h  : boxplot 返回的图元句柄结构体

    % ---------- 参数 ----------
    p = inputParser;
    addParameter(p,'Parent',[],@(x) isempty(x) || isgraphics(x,'axes'));
    addParameter(p,'XLabels',[]);
    addParameter(p,'YLabel','Value',@(x)ischar(x)||isstring(x));
    addParameter(p,'Title','Boxplot by window',@(x)ischar(x)||isstring(x));
    addParameter(p,'Colors',[],@(x)isnumeric(x) && size(x,2)==3);
    addParameter(p,'BoxFaceAlpha',0,@(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    addParameter(p,'Symbol','k+',@(x)ischar(x)||isstring(x));
    addParameter(p,'YScaleFactor',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
    addParameter(p,'YLim',[],@(x) isempty(x) || (isnumeric(x)&&numel(x)==2&&x(1)<x(2)));
    addParameter(p,'AutoTight',true,@(x)islogical(x)&&isscalar(x));
    parse(p, varargin{:});

    ax          = p.Results.Parent;           if isempty(ax), ax = gca; end
    XLabels     = p.Results.XLabels;
    YLabelStr   = char(p.Results.YLabel);
    TitleStr    = char(p.Results.Title);
    Colors      = p.Results.Colors;
    boxAlpha    = p.Results.BoxFaceAlpha;
    symb        = char(p.Results.Symbol);
    scale       = p.Results.YScaleFactor;
    yLimUser    = p.Results.YLim;
    autoTight   = p.Results.AutoTight;

    % ---------- 数据整理 ----------
    [nTrials, nWins] = size(M);
    if isempty(XLabels)
        XLabels = arrayfun(@(k) sprintf('W%d',k), 1:nWins, 'uni', 0);
    else
        assert(numel(XLabels)==nWins, 'XLabels 长度必须等于列数（窗口数）。');
    end

    Yraw = M(:);
    G    = repelem(1:nWins, nTrials).';
    valid = ~isnan(Yraw);
    Yraw = Yraw(valid);
    G    = G(valid);

    % 仅用于显示的缩放（不改变数据本身）
    Y = Yraw * scale;

    % 用带标签 categorical，保持组序
    G = categorical(G, 1:nWins, XLabels);

    % ---------- 颜色准备 ----------
    useColors = ~isempty(Colors);
    if useColors
        if size(Colors,1)==1
            Colors = repmat(Colors, nWins, 1);
        elseif size(Colors,1)~=nWins
            Colors = Colors(mod(0:nWins-1,size(Colors,1))+1, :);
        end
        colors_for_boxplot = Colors;
    else
        colors_for_boxplot = [];
    end

    % ---------- 绘图 ----------
    axes(ax); hold(ax,'on');
    if isempty(colors_for_boxplot)
        h = boxplot(Y, G, 'Symbol', symb, 'Whisker', 1.5);
    else
        h = boxplot(Y, G, 'Symbol', symb, 'Whisker', 1.5, ...
            'Colors', colors_for_boxplot);
    end
    grid(ax,'on');
    title(ax, TitleStr);
    set(ax, 'FontSize', 11, 'XTickLabel', XLabels);

    % y 轴标签：根据 scale 自动附上 ×10^{k}
    ylab = YLabelStr;
    if ~isempty(scale) && isfinite(scale) && scale~=1
        p10 = floor(log10(scale));
        % 纯 10 的幂次（显示更简洁）
        if abs(scale - 10^p10) < 1e-12
            ylab = sprintf('%s (\\times10^{%d})', ylab, p10);
        else
            ylab = sprintf('%s (\\times%.3g)', ylab, scale);
        end
    end
    ylabel(ax, ylab);

    % ---------- 半透明填充（可选） ----------
    if useColors && boxAlpha > 0
        boxes = findobj(ax, 'Tag', 'Box');
        % boxes 的顺序与组别相反，按 X 中心排序为 1..nWins
        cx = arrayfun(@(b) mean(get(b,'XData')), boxes);
        [~, idx] = sort(cx, 'ascend');
        boxes = boxes(idx);

        for k = 1:min(nWins, numel(boxes))
            b = boxes(k);
            x = get(b,'XData'); y = get(b,'YData');
            c = Colors(k,:);
            patch('XData', x, 'YData', y, ...
                  'FaceColor', c, 'FaceAlpha', boxAlpha, ...
                  'EdgeColor', 'none', ...
                  'Parent', ax);
            uistack(b,'top');  % 让箱体边与中位数线保持可见
        end
        meds = findobj(ax,'Tag','Median');
        set(meds,'LineWidth',1.8);
    end

    % ---------- y 轴范围 ----------
    if ~isempty(yLimUser)
        % 用户显式给定（注意：这里是显示后的尺度范围）
        ylim(ax, yLimUser);
    elseif autoTight
        if isempty(Y)
            ylim(ax, [0 1]); % 空数据兜底
        else
            yMin = min(Y); yMax = max(Y);
            if yMax==yMin, yMax = yMin + 1; end
            pad  = 0.08*(yMax - yMin + eps);
            ylim(ax, [yMin-pad, yMax+pad]);
        end
    else
        % 保留当前 y 轴范围
    end
end



% function [ax, h] = boxplot_diy(M, varargin)
% % 对 “试次 × 窗口” 矩阵 M 画箱线图（兼容 subplot）
% %
% % 主要特性：
% % - 支持在指定 axes 上绘图：'Parent', ax
% % - 支持字符串标题/坐标轴：'Title', 'YLabel'（string/char均可）
% % - 自定义颜色：'Colors'（1x3 或 nWins x 3，RGB 0-1）
% % - 半透明箱体：'BoxFaceAlpha'（0~1），需配合 Colors 使用
% % - 其它转发：'Symbol'（离群点标记）
% %
% % 返回：
% %   ax : 使用的坐标轴句柄
% %   h  : boxplot 返回的图元句柄结构体
% 
%     % ---------- 参数 ----------
%     p = inputParser;
%     addParameter(p,'Parent',[],@(x) isempty(x) || isgraphics(x,'axes'));
%     addParameter(p,'XLabels',[]);
%     addParameter(p,'YLabel','Value',@(x)ischar(x)||isstring(x));
%     addParameter(p,'Title','Boxplot by window',@(x)ischar(x)||isstring(x));
%     addParameter(p,'Colors',[],@(x)isnumeric(x) && size(x,2)==3);
%     addParameter(p,'BoxFaceAlpha',0,@(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
%     addParameter(p,'Symbol','k+',@(x)ischar(x)||isstring(x)); % 离群点标记
%     parse(p, varargin{:});
% 
%     ax          = p.Results.Parent;
%     XLabels     = p.Results.XLabels;
%     YLabelStr   = char(p.Results.YLabel);
%     TitleStr    = char(p.Results.Title);
%     Colors      = p.Results.Colors;
%     boxAlpha    = p.Results.BoxFaceAlpha;
%     symb        = char(p.Results.Symbol);
% 
%     if isempty(ax), ax = gca; end
%     axes(ax); 
%     hold(ax,'on');
% 
%     % ---------- 数据整理 ----------
%     [nTrials, nWins] = size(M);
%     Y = M(:);
%     G = repelem(1:nWins, nTrials).'; % 组：窗口编号
% 
%     % 去 NaN
%     valid = ~isnan(Y);
%     Y = Y(valid);
%     G = G(valid);
% 
%     % 分组类别
%     if isempty(XLabels)
%         XLabels = arrayfun(@(k) sprintf('W%d',k), 1:nWins, 'uni', 0);
%     end
%     % 用带标签的 categorical，避免顺序错乱
%     G = categorical(G, 1:nWins, XLabels);
% 
%     % ---------- 颜色准备 ----------
%     useColors = ~isempty(Colors);
%     if useColors
%         if size(Colors,1)==1
%             % 单色扩展为每组一色
%             Colors = repmat(Colors, nWins, 1);
%         elseif size(Colors,1)~=nWins
%             % 长度不匹配则循环使用
%             Colors = Colors(mod(0:nWins-1,size(Colors,1))+1, :);
%         end
%         % boxplot 的 'Colors' 需要 (k x 3)，会循环赋予每个箱线的边颜色
%         colors_for_boxplot = Colors;
%     else
%         colors_for_boxplot = []; % 让 boxplot 用默认色
%     end
% 
%     % ---------- 绘图（画到 ax 上） ----------
%     axes(ax); 
%     hold(ax,'on');
%     
%     if isempty(colors_for_boxplot)
%         h = boxplot(Y, G, 'Symbol', symb, 'Whisker', 1.5);
%     else
%         h = boxplot(Y, G, 'Symbol', symb, 'Whisker', 1.5, ...
%             'Colors', colors_for_boxplot);
%     end
%     grid(ax,'on');
%     ylabel(ax, YLabelStr);
%     title(ax, TitleStr);
%     set(ax, 'FontSize', 11);
%     set(ax, 'XTickLabel', XLabels);
% 
% 
% 
%     % ---------- 半透明填充（可选） ----------
%     % boxplot 默认不填充颜色；如需淡色背景，给每个箱体打补丁。
%     if useColors && boxAlpha > 0
%         % 找到每个箱体的线条对象（Tag='Box'），从其 X/Y 顶点生成 patch
%         boxes = findobj(ax, 'Tag', 'Box');
%         % boxes 返回顺序与组别相反，按位置重排
%         % 通过其 XData 的中位数来排序对应到 1..nWins
%         cx = arrayfun(@(b) mean(get(b,'XData')), boxes);
%         [~, idx] = sort(cx, 'ascend');
%         boxes = boxes(idx);
% 
%         for k = 1:min(nWins, numel(boxes))
%             b = boxes(k);
%             x = get(b,'XData'); y = get(b,'YData');
%             c = Colors(k,:);
%             patch('XData', x, 'YData', y, ...
%                   'FaceColor', c, 'FaceAlpha', boxAlpha, ...
%                   'EdgeColor', 'none', ...
%                   'Parent', ax);
%             % 为了让中位数与边框仍可见，把箱线提到顶层
%             uistack(b,'top');
%         end
% 
%         % 中位数线加粗一点更清晰
%         meds = findobj(ax,'Tag','Median');
%         set(meds,'LineWidth',1.8);
%     end
% end


