function diag = oir_consistency_diagnostics(results, tol)
if nargin<2, tol = 0.05; end   % 允许 5% 相对误差
diag = struct();
for ic = 1:numel(results.cond)
    R = results.cond(ic);
    comps = {'dO12','dO1_2','dO1o2','dO2_1'};
    for k = 1:numel(comps)
        comp = comps{k};
        yspec = R.([comp 'f']);     % 频谱
        tint  = trapz(R.f, yspec) * R.scale;  % 频域积分（带同一scale）
        tval  = R.(comp);            % 时域值
        rel   = abs(tint - tval) / max(1e-9, abs(tval));
        diag(ic).(comp).int    = tint;
        diag(ic).(comp).time   = tval;
        diag(ic).(comp).relerr = rel;
        diag(ic).(comp).pass   = rel <= tol;
        
    end
end
end

