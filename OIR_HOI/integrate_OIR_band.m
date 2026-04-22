function O = integrate_OIR_band(OIR_results_human, band)
% 将频域密度 dO·f 按给定频段积分为 O-information（试次×窗口）
%
% 输入
%   OIR_results_human: 具有 result(file).win(i).O1_2 / O2_1 / O1o2 / O12 的结构
%   band: [] 或 [fmin fmax]
%         - 为空时：全频带积分
%         - 为 [fmin fmax]: 仅在该频段内积分（需要 OIR_results_human.freq 或 .f）
%
% 输出（结构体 O）
%   O.O12, O.O1_2, O.O2_1, O.O1o2, O.residual: 试次×窗口 的矩阵
%   O.freq: 使用的频率轴（若无则为空）
%   O.band: 实际使用的频段（若无则为空）
%
% 备注
% - 若各试次窗口数量不同，自动用 NaN 进行右侧填充
% - residual = O12 - (O1_2 + O2_1 + O1o2)（应接近 0）
% - 若没有频率轴且指定了 band，将按频点索引近似截取（给出警告）

    % 频率轴
    if isfield(OIR_results_human, 'freq')
        f = OIR_results_human.freq(:);
    elseif isfield(OIR_results_human, 'f')
        f = OIR_results_human.f(:);
        disp(size(f))
    else
        f = [];
    end

    use_full_band = isempty(band);

    % 定位需要积分的频点索引
    if ~use_full_band
        if ~isempty(f)
            fmin = band(1); fmax = band(2);
            if fmin > fmax, error('band 必须满足 [fmin fmax] 且 fmin<=fmax'); end
            idx = find(f >= fmin & f <= fmax);
            if isempty(idx), warning('指定频段内无频点；结果将全为 NaN。'); end
        else
            % 没有频率轴：退化为按索引范围近似
            warning('未提供频率轴(freq/f)，band 将按频点索引近似截取。');
            idx = []; % 稍后按长度构造
        end
    else
        idx = []; % 全频带时稍后忽略
    end

    % 尺度
    assert(isfield(OIR_results_human, 'result'), '缺少 result 字段。');
    nFiles = numel(OIR_results_human.result);
    nWins_each = arrayfun(@(x) numel(x.win), OIR_results_human.result);
    maxWins = max(nWins_each);

    % 预分配（试次×窗口）
    O12   = nan(nFiles, maxWins);
    O1_2  = nan(nFiles, maxWins);
    O2_1  = nan(nFiles, maxWins);
    O1o2  = nan(nFiles, maxWins);
    RES   = nan(nFiles, maxWins);

    for fi = 1:nFiles
        wins = OIR_results_human.result(fi).win;
        for wi = 1:numel(wins)
            d12  = vec(wins(wi).O12);
            d1_2 = vec(wins(wi).O1_2);
            d2_1 = vec(wins(wi).O2_1);
            dco  = vec(wins(wi).O1o2);

            L = unique([numel(d12), numel(d1_2), numel(d2_1), numel(dco)]);
            if numel(L)~=1
                error('谱长不一致：file %d win %d', fi, wi);
            end
            L = L(1);

            % 确定积分索引
            if use_full_band
                useIdx = 1:L;
            else
                if isempty(f)
                    % 仅有长度，按“索引对应频点”近似（fmin..fmax 映射为相对比例）
                    useIdx = map_band_to_index(band, L);
                else
                    useIdx = idx;
                end
            end
            if isempty(useIdx), I12=NaN; I1_2=NaN; I2_1=NaN; Io=NaN;
            else
                if isempty(f) || ~use_full_band && isempty(f)
                    % 按索引积分（相对尺度）
                    x = (1:L).';
                    I12  = trapz(x(useIdx),  d12(useIdx));
                    I1_2 = trapz(x(useIdx), d1_2(useIdx));
                    I2_1 = trapz(x(useIdx), d2_1(useIdx));
                    Io   = trapz(x(useIdx),  dco(useIdx));
                else
                    % 按真实频率积分（物理 Hz 尺度）
                    I12  = trapz(f(useIdx),  d12(useIdx));
                    I1_2 = trapz(f(useIdx), d1_2(useIdx));
                    I2_1 = trapz(f(useIdx), d2_1(useIdx));
                    Io   = trapz(f(useIdx),  dco(useIdx));
                end
            end

            O12(fi,wi)  = I12;
            O1_2(fi,wi) = I1_2;
            O2_1(fi,wi) = I2_1;
            O1o2(fi,wi) = Io;
            RES(fi,wi)  = I12 - (I1_2 + I2_1 + Io);
        end
    end

    % 输出
    O = struct();
    O.O12      = O12;
    O.O1_2     = O1_2;
    O.O2_1     = O2_1;
    O.O1o2     = O1o2;
    O.residual = RES;
    O.freq     = f;
    O.band     = ifelse(use_full_band, [], band);
end

% ------- 小工具 --------
function y = vec(x)
    if isrow(x), y = x.'; else, y = x; end
end

function idx = map_band_to_index(band, L)
% 在没有频率轴时，用索引近似：把 band 的相对位置线性映射到 [1..L]
% 假设原频带是 [0..1] 的相对坐标（仅用于粗略子带比较）
    fmin = max(0, min(1, band(1)));
    fmax = max(0, min(1, band(2)));
    if fmin > fmax, tmp=fmin; fmin=fmax; fmax=tmp; end
    i1 = max(1, round(1 + fmin*(L-1)));
    i2 = max(1, round(1 + fmax*(L-1)));
    idx = i1:i2;
end

function v = ifelse(c,a,b)
    if c, v=a; else, v=b; end
end

