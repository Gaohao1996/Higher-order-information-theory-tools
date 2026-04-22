function draw_group_contours(data2d, group_id, n_groups)

cols = lines(n_groups);

xg = linspace(min(data2d(:,1))-0.5, max(data2d(:,1))+0.5, 120);
yg = linspace(min(data2d(:,2))-0.5, max(data2d(:,2))+0.5, 120);
[Xg, Yg] = meshgrid(xg, yg);
XY = [Xg(:), Yg(:)];

for g = 1:n_groups
    idx = group_id == g;
    D = data2d(idx, :);

    if size(D,1) < 5
        continue;
    end

    mu = mean(D,1);
    S = cov(D);

    if rank(S) < 2
        S = S + 1e-6 * eye(2);
    end

    try
        Zpdf = mvnpdf(XY, mu, S);
        Zpdf = reshape(Zpdf, size(Xg));

        lev = max(Zpdf(:)) * 0.35;
        contour(Xg, Yg, Zpdf, [lev lev], ...
            'Color', cols(g,:), 'LineWidth', 1.8);
    catch
        % 若 mvnpdf 因数值问题失败，就跳过
    end
end

end