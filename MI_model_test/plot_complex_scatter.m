function plot_complex_scatter(x, y, Z, group_id, color_mode)

switch lower(color_mode)
    case 'continuous'
        scatter(x, y, 14, Z, 'filled', ...
            'MarkerFaceAlpha', 0.70, 'MarkerEdgeAlpha', 0.70);
        colormap(parula);
        colorbar;

    case 'group'
        ng = max(group_id);
        cols = lines(ng);
        for g = 1:ng
            idx = group_id == g;
            scatter(x(idx), y(idx), 16, cols(g,:), 'filled', ...
                'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
            hold on;
        end
    otherwise
        error('未知 color_mode');
end

end