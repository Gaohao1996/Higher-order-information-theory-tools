function plot_density2d(x, y)

nGrid = 60;
xv = linspace(min(x), max(x), nGrid);
yv = linspace(min(y), max(y), nGrid);
[Xg, Yg] = meshgrid(xv, yv);

xy = [x(:), y(:)];
[f, ~] = ksdensity(xy, [Xg(:), Yg(:)]);
Fg = reshape(f, size(Xg));

contour(Xg, Yg, Fg, 8, 'LineWidth', 1.2);
end