function h = plot_cov_ellipse(X, color_use, nsig, line_style, line_width)
% Plot covariance ellipse for 2D data
%
% X         : N x 2 data
% color_use : RGB color
% nsig      : ellipse radius in std units, e.g. 1, 2
% line_style: '-', '--', ':'
% line_width: scalar

    if nargin < 3 || isempty(nsig)
        nsig = 2;
    end
    if nargin < 4 || isempty(line_style)
        line_style = '-';
    end
    if nargin < 5 || isempty(line_width)
        line_width = 1.5;
    end

    mu = mean(X, 1);
    C = cov(X);

    if any(isnan(C(:))) || rank(C) < 2
        h = gobjects(1);
        return;
    end

    [V, D] = eig(C);
    t = linspace(0, 2*pi, 200);
    circle = [cos(t); sin(t)];

    % ellipse in data space
    ellipse = (V * sqrt(D) * nsig * circle)';
    ellipse = ellipse + mu;

    h = plot(ellipse(:,1), ellipse(:,2), ...
        'Color', color_use, ...
        'LineStyle', line_style, ...
        'LineWidth', line_width);
end