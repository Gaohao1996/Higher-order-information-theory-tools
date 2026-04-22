function O = integrate_OIR_band_new(results, f, band, varargin)
% 将 OIR 频谱密度按给定频段积分为 O-information（试次×窗口）
%
% 输入
%   results : 结构数组 = OIR_results_human.results
%             每个元素含 .win(w).O1_2 / .O2_1 / .O1o2 / .O12 ，均为列或行向量
%   f       : 频率轴 (Hz)，长度与各 O* 一致，单调递增
%   band    : [] 或 [fmin fmax]（左闭右开：f >= fmin & f < fmax）
%
% 可选参数（Name-Value）：
%   'BandwidthNormalize' (false)  true 则除以带宽，得到带内平均强度
%   'Fmin' ([])                   设为 0.5 或 1 可排除 DC/超低频影响
%
% 输出
%   O: 结构体，含 O12 / O1_2 / O2_1 / O1o2 / residual （试次×窗口）
%      以及 O.freq (=f) 与 O.band (=band)

    p = inputParser;
    addParameter(p,'BandwidthNormalize',false,@islogical);
    addParameter(p,'Fmin',[],@(x) isempty(x) || isscalar(x));
    parse(p,varargin{:});
    bwNorm  = p.Results.BandwidthNormalize;
    fminCut = p.Results.Fmin;

    f = f(:);
    assert(isvector(f) && issorted(f),'f 必须为单调递增向量');
    nFiles   = numel(results);
    nWinsEach = arrayfun(@(s) numel(s.win), results);
    maxWins  = max(nWinsEach);

    % 频段掩码：左闭右开，避免端点重复
    if isempty(band)
        idx = 1:numel(f);
        eff_bandwidth = f(end) - f(1);
    else
        f_lo = band(1); f_hi = band(2);
        assert(f_lo <= f_hi, 'band 必须满足 [fmin fmax] 且 fmin<=fmax');
        if ~isempty(fminCut)
            f_lo = max(f_lo, fminCut);
        end
        idx = find(f >= f_lo & f < f_hi);
        if isempty(idx)
            warning('指定频段 [%g, %g) Hz 无频点。', f_lo, f_hi);
        end
        if ~isempty(idx)
            eff_bandwidth = f(idx(end)) - f(idx(1));
        else
            eff_bandwidth = NaN;
        end
    end

    % 预分配（试次×窗口）
    O12  = nan(nFiles, maxWins);
    O1_2 = nan(nFiles, maxWins);
    O2_1 = nan(nFiles, maxWins);
    O1o2 = nan(nFiles, maxWins);
    RES  = nan(nFiles, maxWins);

    for fi = 1:nFiles
        for wi = 1:numel(results(fi).win)
            d12  = vec(results(fi).win(wi).O12);
            d1_2 = vec(results(fi).win(wi).O1_2);
            d2_1 = vec(results(fi).win(wi).O2_1);
            dco  = vec(results(fi).win(wi).O1o2);

            % 基于 Hz 的 trapz 积分
            if isempty(idx)
                I12=NaN; I1_2=NaN; I2_1=NaN; Io=NaN;
            else
                I12  = trapz(f(idx),  d12(idx));
                I1_2 = trapz(f(idx), d1_2(idx));
                I2_1 = trapz(f(idx), d2_1(idx));
                Io   = trapz(f(idx),  dco(idx));

                if bwNorm && ~isnan(eff_bandwidth) && eff_bandwidth>0
                    I12  = I12  / eff_bandwidth;
                    I1_2 = I1_2 / eff_bandwidth;
                    I2_1 = I2_1 / eff_bandwidth;
                    Io   = Io   / eff_bandwidth;
                end
            end

            O12(fi,wi)  = I12;
            O1_2(fi,wi) = I1_2;
            O2_1(fi,wi) = I2_1;
            O1o2(fi,wi) = Io;
            RES(fi,wi)  = I12 - (I1_2 + I2_1 + Io);
        end
    end

    O = struct();
    O.O12      = O12;
    O.O1_2     = O1_2;
    O.O2_1     = O2_1;
    O.O1o2     = O1o2;
    O.residual = RES;
    O.freq     = f;
    O.band     = band;
end

% helper
function y = vec(x), if isrow(x), y = x.'; else, y = x; end, end


