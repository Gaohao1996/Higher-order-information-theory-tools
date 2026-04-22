function results = oir_flexible(data_cells, fs, opts)
% Flexible OIR analysis across conditions with arbitrary block grouping.
% ---------------------------------------------------------------------
% data_cells : 1xC cell, each is [N x P] (rows=samples, cols=variables)
% fs         : sampling rate in Hz
% opts       : struct with fields:
%   .var_order  : 1xP' column order applied to each condition (default 1:P)
%   .Mv         : 1xB block sizes; sum(Mv)==P'  (REQUIRED)
%   .iM         : vector of BLOCK indices participating (default 1:B)  <-- 块索引！
%   .j          : scalar TARGET BLOCK index (1..B, j ∈ iM)            <-- 块索引！
%   .nfft       : #points for [0, fs/2] spectrum (default 1000)
%   .pmax       : max VAR order for AIC selection (default 20); ignored if .order provided
%   .order      : fixed VAR order (optional)
%   .bands      : Kx2 cell, e.g. {'alpha',[8 12]; 'beta',[13 30]} (Hz)
%   .peakBand   : struct with .name, .chan, .width, .search  (per-condition peak-centered band)
%                 .chan is a COLUMN index AFTER var_order (not block index)
%   .plot       : true/false (default true)
%   .cond_names : 1xC cellstr (default {'Cond1', ...})
%
% RETURNS:
% results.info      : copy of opts + notes
% results.cond(k)   : struct with fields:
%   f, dO12f, dO1_2f, dO1o2f, dO2_1f, dO12, dO1_2, dO1o2, dO2_1, scale, bands

% ---------- sanity checks / defaults ----------
C = numel(data_cells);
assert(C>=1, 'data_cells must be a non-empty cell array');

P0 = size(data_cells{1},2);
for ic=2:C
    assert(size(data_cells{ic},2)==P0, 'All conditions must have same #columns');
end

if ~isfield(opts,'var_order') || isempty(opts.var_order)
    opts.var_order = 1:P0;
end
Pp = numel(opts.var_order);              % #columns after reorder

if ~isfield(opts,'Mv') || isempty(opts.Mv)
    error('opts.Mv (block sizes) is required');
end
assert(sum(opts.Mv)==Pp, 'sum(Mv) must equal number of variables after var_order');

B = numel(opts.Mv);                      % #blocks
if ~isfield(opts,'iM') || isempty(opts.iM)
    opts.iM = 1:B;                       % iM 是块索引
else
    assert(all(ismember(opts.iM,1:B)), 'iM must be BLOCK indices in 1..B');
end
if ~isfield(opts,'j') || isempty(opts.j)
    error('opts.j (target BLOCK index) is required');
end
assert(ismember(opts.j, 1:B), 'j must be a BLOCK index in 1..B');
assert(ismember(opts.j, opts.iM), 'j must be one of the participating blocks (iM)');

if ~isfield(opts,'nfft')  || isempty(opts.nfft),  opts.nfft = 1000; end
if ~isfield(opts,'pmax')  || isempty(opts.pmax),  opts.pmax = 20;   end
if ~isfield(opts,'plot')  || isempty(opts.plot),  opts.plot = true; end
if ~isfield(opts,'cond_names') || isempty(opts.cond_names)
    opts.cond_names = arrayfun(@(k) sprintf('Cond%d',k), 1:C, 'uni',0);
end

fixedBands = {};
if isfield(opts,'bands') && ~isempty(opts.bands)
    fixedBands = opts.bands;
    assert(iscell(fixedBands) && size(fixedBands,2)==2, ...
        'opts.bands must be Kx2 cell: {name,[f1 f2]; ...}');
end

usePeakBand = isfield(opts,'peakBand') && ~isempty(opts.peakBand);
if usePeakBand
    pb = opts.peakBand;
    need = {'name','chan','width','search'};
    for t=1:numel(need), assert(isfield(pb,need{t}), 'peakBand.%s required', need{t}); end
    assert(pb.chan>=1 && pb.chan<=Pp, 'peakBand.chan must be a column index AFTER var_order (1..P'')');
end

% ---------- preallocate results ----------
results = struct();
results.info = opts;
results.info.note = 'IMPORTANT: iM and j are BLOCK indices (1..numel(Mv)).';

cond_template = struct( ...
    'f',[], ...
    'dO12f',[], 'dO1_2f',[], 'dO1o2f',[], 'dO2_1f',[], ...
    'dO12',nan, 'dO1_2',nan, 'dO1o2',nan, 'dO2_1',nan, ...
    'scale',1, ...
    'bands',struct() );
results.cond = repmat(cond_template, C, 1);

% ---------- main loop ----------
for ic = 1:C
    S0 = data_cells{ic};
    S  = S0(:, opts.var_order);                % reorder columns
    S  = S - mean(S,1,'omitnan');              % remove mean
    
    % VAR order
    if isfield(opts,'order') && ~isempty(opts.order)
        order = opts.order;
    else
        order = oir_mosVAR(S', opts.pmax);     % AIC selection
    end
    
    %test the influence by order change
    %  order = 8;
%     disp(order)
    
    % Fit VAR
    [Am,Su] = oir_idVAR(S', order);
    % VAR -> ISS
    [A,Ciss,K,~] = oir_ar2iss(Am);
    
    % Compute OIR (block indices!)

    %% print drivers and target(options)
    % fprintf('OIR group ij = [%s], target j = %d, drivers = [%s]\n', ...
    % num2str(opts.iM), opts.j, num2str(setdiff(opts.iM, opts.j)));

    out = oir_deltaO(A,Ciss,K,Su,opts.Mv,opts.iM,opts.j,fs,opts.nfft);
    
    % Collect spectra & time scalars
    R = struct();
    R.f       = out.freq(:);
    R.dO12f   = out.dO12f(:);
    R.dO1_2f  = out.dO1_2f(:);
    R.dO1o2f  = out.dO1o2f(:);
    R.dO2_1f  = out.dO2_1f(:);
    R.dO12    = out.dO12;
    R.dO1_2   = out.dO1_2;
    R.dO1o2   = out.dO1o2;
    R.dO2_1   = out.dO2_1;
    
    
    %     % Scale so that integral of spectrum matches time-domain scalar
    %     tot_spec = trapz(R.f, R.dO12f);
    %     if ~isfinite(tot_spec) || abs(tot_spec) < eps
    %         R.scale = 1;
    %     else
    %         R.scale = R.dO12 / tot_spec;
    %     end
    
    % ---- print ----
    I_tot   = trapz(R.f, R.dO12f);
    ratio_t = I_tot / (R.dO12 + eps);
    
    I_dir   = trapz(R.f, R.dO1_2f) / (R.dO1_2 + eps);
    I_cpl   = trapz(R.f, R.dO1o2f) / (R.dO1o2 + eps);
    I_rev   = trapz(R.f, R.dO2_1f) / (R.dO2_1 + eps);
    
    % fprintf('[OIR check | %s] total=%.3f dir=%.3f cpl=%.3f rev=%.3f\n', ...
    %     opts.cond_names{ic}, ratio_t, I_dir, I_cpl, I_rev);
    
    if any(~isfinite([ratio_t I_dir I_cpl I_rev])) || ...
            any(abs([ratio_t I_dir I_cpl I_rev]-1) > 0.03)
        warning('OIR:UnitMismatch', ...
            '频域积分与时域标量不够一致，请检查上游频谱规范化。');
    end
    
    R.scale = 1;  % 不再使用外部 scale
    
    
    % ---- band integrations ----
    bandsResult = struct();
    
    % (1) fixed bands
    for kb = 1:size(fixedBands,1)
        bname  = fixedBands{kb,1};
        brange = fixedBands{kb,2};
        mask   = (R.f >= brange(1)) & (R.f <= brange(2));
        %         bandsResult.(bname).range = brange;
        %         bandsResult.(bname).dO12  = R.scale*trapz(R.f(mask), R.dO12f(mask));
        %         bandsResult.(bname).dO1_2 = R.scale*trapz(R.f(mask), R.dO1_2f(mask));
        %         bandsResult.(bname).dO1o2 = R.scale*trapz(R.f(mask), R.dO1o2f(mask));
        %         bandsResult.(bname).dO2_1 = R.scale*trapz(R.f(mask), R.dO2_1f(mask));
        bandsResult.(bname).dO12  = trapz(R.f(mask), R.dO12f(mask));
        bandsResult.(bname).dO1_2 = trapz(R.f(mask), R.dO1_2f(mask));
        bandsResult.(bname).dO1o2 = trapz(R.f(mask), R.dO1o2f(mask));
        bandsResult.(bname).dO2_1 = trapz(R.f(mask), R.dO2_1f(mask));
        
    end
    
    % (2) peak-centered band per condition
    if usePeakBand
        pb = opts.peakBand;
        x  = S(:, pb.chan);                 % column after var_order
        % Welch PSD for peak search
        nseg = min(length(x), 4096); nseg = 2^floor(log2(nseg));
        if nseg < 256, nseg = 256; end
        [pxx,ff] = pwelch(x, hamming(nseg), round(0.5*nseg), max(2048,nseg), fs);
        msk = (ff >= pb.search(1)) & (ff <= pb.search(2));
        if any(msk)
            [~, idx] = max(pxx(msk));
            fsub  = ff(msk);
            fpeak = fsub(idx);
        else
            fpeak = mean(pb.search);        % fallback
        end
        brange = [max(0, fpeak - pb.width), min(fs/2, fpeak + pb.width)];
        mask   = (R.f >= brange(1)) & (R.f <= brange(2));
        
        bandsResult.(pb.name).range = brange;
        bandsResult.(pb.name).fpeak = fpeak;
        bandsResult.(pb.name).dO12  = R.scale*trapz(R.f(mask), R.dO12f(mask));
        bandsResult.(pb.name).dO1_2 = R.scale*trapz(R.f(mask), R.dO1_2f(mask));
        bandsResult.(pb.name).dO1o2 = R.scale*trapz(R.f(mask), R.dO1o2f(mask));
        bandsResult.(pb.name).dO2_1 = R.scale*trapz(R.f(mask), R.dO2_1f(mask));
    end
    
    R.bands = bandsResult;
    
    % safe struct assignment (avoid "不同结构体之间下标赋值" error)
    results.cond(ic).f       = R.f;
    results.cond(ic).dO12f   = R.dO12f;
    results.cond(ic).dO1_2f  = R.dO1_2f;
    results.cond(ic).dO1o2f  = R.dO1o2f;
    results.cond(ic).dO2_1f  = R.dO2_1f;
    results.cond(ic).dO12    = R.dO12;
    results.cond(ic).dO1_2   = R.dO1_2;
    results.cond(ic).dO1o2   = R.dO1o2;
    results.cond(ic).dO2_1   = R.dO2_1;
    results.cond(ic).scale   = R.scale;
    results.cond(ic).bands   = R.bands;
end

% ---------- optional plotting ----------
% if opts.plot
%     K = C;
%     figure('Color','w','Name','OIR Spectral Decomposition');
%     tl = tiledlayout(3, K, 'TileSpacing','compact','Padding','compact');
%     
%     ymin = +inf; ymax = -inf;
%     for ic=1:C
%         R = results.cond(ic);
%         ymin = min([ymin; R.dO12f(:); R.dO1_2f(:); R.dO1o2f(:); R.dO2_1f(:)]);
%         ymax = max([ymax; R.dO12f(:); R.dO1_2f(:); R.dO1o2f(:); R.dO2_1f(:)]);
%     end
%     pad = 0.06*(ymax - ymin + eps);
%     
%     for ic=1:C
%         R = results.cond(ic);
%         
%         % Row 1: total
%         nexttile(ic);
%         plot(R.f, R.dO12f, 'LineWidth', 1.2); grid on; hold on;
%         title(opts.cond_names{ic});
%         ylabel('\delta_{total}(f)');
%         xlim([0 fs/2]); ylim([ymin-pad, ymax+pad]); yline(0,'k:');
%         shade_bands_if_any(R);
%         
%         % Row 2: dir + coupling
%         nexttile(K+ic);
%         plot(R.f, R.dO1_2f, 'LineWidth', 1.2); hold on;
%         plot(R.f, R.dO1o2f, 'LineWidth', 1.2); grid on;
%         if ic==1, ylabel('dir→j  &  coupling'); end
%         legend({'dir (sources→j)','coupling'}, 'Location','best'); legend boxoff;
%         xlim([0 fs/2]); ylim([ymin-pad, ymax+pad]); yline(0,'k:');
%         shade_bands_if_any(R);
%         
%         % Row 3: reverse j→sources
%         nexttile(2*K+ic);
%         plot(R.f, R.dO2_1f, 'LineWidth', 1.2); grid on;
%         if ic==1, ylabel('reverse (j→sources)'); end
%         xlabel('Frequency (Hz)'); xlim([0 fs/2]); ylim([ymin-pad, ymax+pad]); yline(0,'k:');
%         shade_bands_if_any(R);
%     end
% end
end

% ---- helper to shade bands on current axes ----
function shade_bands_if_any(R)
if ~isfield(R,'bands') || isempty(fieldnames(R.bands)), return; end
yl = ylim;
bn = fieldnames(R.bands);
for k = 1:numel(bn)
    rg = R.bands.(bn{k}).range;
    patch([rg(1) rg(2) rg(2) rg(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0 0 0], 'FaceAlpha',0.06, 'EdgeColor','none');
    text(mean(rg), yl(2)-0.04*(yl(2)-yl(1)), bn{k}, ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end
uistack(findobj(gca,'Type','line'),'top');
end
