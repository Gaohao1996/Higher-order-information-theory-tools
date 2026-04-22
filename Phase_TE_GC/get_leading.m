% 用局部函数体实现（兼容老版 MATLAB 无匿名嵌套）：
function pos_keep = get_leading(sig_idx, total_len)
    % total_len = numel(drivers_*), 位置 1=i, 2=j，真正新增从 3 开始
    if total_len <= 2
        pos_keep = [];
        return;
    end
    pos_all = 3:total_len;                          % 只看新增通道的位置序列
    sig_flags = ismember(pos_all, sig_idx(:).');    % 哪些位置被判为显著
    k = find(~sig_flags, 1, 'first');               % 第一处“非显著”
    if isempty(k)
        pos_keep = pos_all;                          % 全部显著
    else
        lead = k - 1;                                % 连续显著前缀长度
        if lead <= 0
            pos_keep = [];                           % 第一个就是非显著 -> 空
        else
            pos_keep = pos_all(1:lead);
        end
    end
end
