function te = CTE_iteration(x_sub, p, tau, i_sub, j_sub, delayZ, phase)
% CTE_iteration
% ------------------------------------------------------------
% 计算条件传递熵 TE(X -> Y | Z)
%
% 输入:
%   x_sub : [T x nSel] 当前用于计算的子矩阵
%   p     : 模型阶数 / embedding order
%   tau   : delay
%   i_sub : x_sub 中 driver 对应的列号，可以是向量
%   j_sub : x_sub 中 target 对应的列号，可以是向量
%   phase : 1 表示相位数据, 0 表示普通连续数据
%
% 说明:
%   - driver = x_sub(:, i_sub)
%   - target = x_sub(:, j_sub)
%   - 剩余列自动作为条件变量 Z
%
% 输出:
%   te    : conditional_TE_delay 返回的 TE 值
% ------------------------------------------------------------
if nargin < 7
    delayZ = [];
end

[~, nSel] = size(x_sub);
all_idx = 1:nSel;

% 转成行向量，避免后续 setdiff / intersect 出现方向问题
i_sub = i_sub(:)';
j_sub = j_sub(:)';

% -------- 输入合法性检查 --------
if isempty(i_sub)
    error('CTE_iteration: i_sub cannot be empty.');
end

if isempty(j_sub)
    error('CTE_iteration: j_sub cannot be empty.');
end

if any(i_sub < 1) || any(i_sub > nSel)
    error('CTE_iteration: i_sub out of range.');
end

if any(j_sub < 1) || any(j_sub > nSel)
    error('CTE_iteration: j_sub out of range.');
end

if ~isempty(intersect(i_sub, j_sub))
    error('CTE_iteration: driver and target indices overlap.');
end

% -------- 提取 X, Y, Z --------
X = x_sub(:, i_sub);
Y = x_sub(:, j_sub);

z_sub = setdiff(all_idx, [i_sub j_sub]);

if isempty(z_sub)
    Z = [];
else
    Z = x_sub(:, z_sub);
end

te = conditional_TE_delay(X, Y, Z, p, tau, delayZ, phase);
end



% function te = CTE_iteration(x, p, tau, i, j, phase)
% % x: [T x nSel] 由 TE_syn_red 传入的子集矩阵 x(:,dt)
% % i,j: 仍然是"原始通道编号"
% % phase: 1 -> 相位数据（弧度），在 conditional_TE_delay 内部展开
% %        0 -> 普通连续数据
% 
%     [~, nSel] = size(x);
% 
%     % --- 关键：把原始索引 i/j 映射到当前子矩阵 x(:,dt) 中的位置 ---
%     % TE_syn_red 传入的是 x(:,dt)，dt 的前两列通常是 [i j]，后面是条件变量
%     % 因此这里最稳妥：默认 driver 在第1列，target 在第2列
%     % 同时做一个保护：如果 nSel<2 直接报错
%     if nSel < 2
%         error('CTE_iteration: x must contain at least [driver target] columns.');
%     end
% 
%     X = x(:, i);
%     Y = x(:, j);
% 
%     if nSel <= 2
%         Z = [];
%     else
%         Z = x(:, 3:end);  % 剩余列全部作为条件变量（保持你原先迭代含义）
%     end
% 
%     deltaZ = 0; % 与你原来一致
%     te = conditional_TE_delay(X, Y, Z, p, tau, deltaZ, phase);
% end






