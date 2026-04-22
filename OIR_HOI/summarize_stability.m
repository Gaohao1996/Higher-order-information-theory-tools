function S = summarize_stability(ReCell, labels)
% ReCell{c}: trials×K（已按秒归一化并重采样）
    nC = numel(ReCell);
    S = struct([]);
    for c = 1:nC
        X = ReCell{c};
        K = size(X,2);
        med = median(X,'omitnan');          % 1×K
        q1  = quantile(X,0.25,1);           % 1×K
        q3  = quantile(X,0.75,1);           % 1×K
        iqr = q3 - q1;                      % 1×K
        % Tukey 异常值率
        lo = q1 - 1.5*iqr;
        hi = q3 + 1.5*iqr;
        isOut = (X < lo) | (X > hi);   
        % 计算每个 bin 的样本数（排除 NaN）
        validN = sum(~isnan(X),1);
        bin_rate = sum(isOut,1) ./ max(validN,1);   % 1×K

        eps = 1e-6;
        cvIQR = iqr./(abs(med)+eps);        % 1×K
        S(c).label = labels{c};
        S(c).median_overall = median(med,'omitnan');
        S(c).IQR_overall    = median(iqr,'omitnan');
        S(c).CV_IQR_score   = median(cvIQR,'omitnan');   % better with small val
        S(c).Outlier_rate   = mean(bin_rate, 'omitnan'); % median(outRate,'omitnan'); % median val version
    end
    % 打印排序
    [~,idx] = sort([S.CV_IQR_score]);
    fprintf('Rank by stability (CV_IQR):\n');
    for r=1:nC
        c = idx(r);
        fprintf('%d) %s | CV_IQR=%.3g | IQR=%.3g | OutRate=%.3g | med=%.3g\n',...
            r,S(c).label,S(c).CV_IQR_score,S(c).IQR_overall,S(c).Outlier_rate,S(c).median_overall);
    end
    
      % ------- 绘图部分 -------
    figure('Name','Stability Summary','Units','normalized','Position',[0.1 0.1 0.8 0.6]);

    % 数据收集
    IQR_vals  = [S.IQR_overall];
    CV_vals   = [S.CV_IQR_score];
    MED_vals  = [S.median_overall];
    OUT_vals  = [S.Outlier_rate];

    x = 1:nC;

    % 子图1：IQR（稳健性）
    subplot(1,4,1);
    bar(x, IQR_vals);
    title('IQR (Stability)');
    set(gca,'XTick',x,'XTickLabel',labels);
    ylabel('IQR');
    grid on;

    % 子图2：CV-IQR
    subplot(1,4,2);
    bar(x, CV_vals);
    title('CV-IQR Score');
    set(gca,'XTick',x,'XTickLabel',labels);
    ylabel('CV IQR');
    grid on;

    % 子图3：Median（效应量）
    subplot(1,4,3);
    bar(x, MED_vals);
    title('Median Effect');
    set(gca,'XTick',x,'XTickLabel',labels);
    ylabel('Median');
    grid on;

    % 子图4：Outlier rate
    subplot(1,4,4);
    bar(x, OUT_vals);
    title('Outlier Rate');
    set(gca,'XTick',x,'XTickLabel',labels);
    ylabel('Rate');
    grid on;

    sgtitle(sprintf('Window Parameter Comparison (Stability & Effect) K=%d',K));
    
end

