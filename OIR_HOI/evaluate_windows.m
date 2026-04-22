function S = evaluate_windows(ReCell, labels)
    nC = numel(ReCell);
    S = struct([]);

    for c = 1:nC
        X = ReCell{c};      % trials×K
        med_k = median(X,1,'omitnan');  
        q1 = quantile(X,0.25,1);
        q3 = quantile(X,0.75,1);
        iqr_k = q3 - q1;

        S(c).label = labels{c};
        S(c).median_effect = median(med_k,'omitnan');   % 效应量（越大越强）
        S(c).stability = median(iqr_k,'omitnan');       % 稳健性（越小越好）
    end

    % 排序打印
    fprintf("\nRanking by stability (IQR median):\n");
    [~,idx1]=sort([S.stability]);
    for r=1:nC
        c=idx1(r);
        fprintf("%d) %s | stability=%.3f | effect=%.3f\n",...
            r,S(c).label,S(c).stability,S(c).median_effect);
    end

    fprintf("\nRanking by effect (median):\n");
    [~,idx2]=sort([S.median_effect],'descend');
    for r=1:nC
        c=idx2(r);
        fprintf("%d) %s | effect=%.3f | stability=%.3f\n",...
            r,S(c).label,S(c).median_effect,S(c).stability);
    end
end


