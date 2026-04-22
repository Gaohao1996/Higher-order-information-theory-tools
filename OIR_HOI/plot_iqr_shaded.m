function h = plot_iqr_shaded(x, data, color, alpha)
% data: [Ntrials x Nfreq]
% x   : [Nfreq x 1] 或 [1 x Nfreq]
if nargin < 3 || isempty(color), color = [0.3 0.5 0.9]; end
if nargin < 4, alpha = 0.3; end

% --- 统一形状 ---
x   = x(:)';                       % 行向量 [1 x M]
if size(data,2) ~= numel(x)
    error('data 的列数(%d)必须等于 x 的长度(%d)。', size(data,2), numel(x));
end

% 分位数 & 中位数（按列）
q25 = squeeze(quantile(data, 0.25, 1)); q25 = q25(:)';   % [1 x M]
q75 = squeeze(quantile(data, 0.75, 1)); q75 = q75(:)';   % [1 x M]
med = squeeze(median(data, 1));        med = med(:)';    % [1 x M]

% --- 阴影 + 中位数 ---
fill([x fliplr(x)], [q75 fliplr(q25)], color, ...
     'FaceAlpha', alpha, 'EdgeColor', 'none'); 
hold on;
h = plot(x, med, 'Color', color, 'LineWidth', 2);
end
