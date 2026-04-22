function out = plot_tf_feature_joint_structure(X, Z, varargin)

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;
addParameter(p, 'mode', 'auto');                 % 'auto' | 'univariate' | 'joint'
addParameter(p, 'nbins', 3);
addParameter(p, 'group_method', 'quantile');     % 'quantile' | 'equal' | 'kmeans'
addParameter(p, 'z_dim_method', 'direct');       % 'direct' | 'pca'
addParameter(p, 'z_dim_idx', 1);
addParameter(p, 'x_method', 'direct');           % 'direct' | 'pca'
addParameter(p, 'x_dims', [1 2]);

% separate control
addParameter(p, 'copula_norm_X', false);
addParameter(p, 'copula_norm_Z', false);

addParameter(p, 'plot_type_1d', 'hist');         % 'hist' | 'box'
addParameter(p, 'plot_type_2d', 'overlay_scatter'); % 'overlay_scatter' | 'subplot_scatter' | 'contour'
addParameter(p, 'title_str', '');
parse(p, varargin{:});
prm = p.Results;

% -------------------------
% Basic checks
% -------------------------
if size(X,1) ~= size(Z,1)
    error('X and Z must have the same number of samples.');
end

X = double(X);
Z = double(Z);
N = size(X,1); %#ok<NASGU>

valid_idx = all(~isnan(X),2) & all(~isnan(Z),2);
X = X(valid_idx,:);
Z = Z(valid_idx,:);

% -------------------------
% Optional copula normalization (separate control for X and Z)
% -------------------------
if prm.copula_norm_X
    Xp = copnorm(X);
    x_space_str = 'copula';
else
    Xp = X;
    x_space_str = 'raw';
end

if prm.copula_norm_Z
    Zp = copnorm(Z);
    z_space_str = 'copula';
else
    Zp = Z;
    z_space_str = 'raw';
end

% -------------------------
% Build scalar Z for grouping
% -------------------------
pcaZ = [];
if size(Zp,2) == 1
    Z_group_base = Zp(:,1);
else
    switch lower(prm.z_dim_method)
        case 'direct'
            Z_group_base = Zp(:, prm.z_dim_idx);
        case 'pca'
            [coeffZ, scoreZ, latentZ, ~, explainedZ] = pca(Zp);
            Z_group_base = scoreZ(:,1);
            pcaZ.coeff = coeffZ;
            pcaZ.score = scoreZ;
            pcaZ.latent = latentZ;
            pcaZ.explained = explainedZ;
        otherwise
            error('Unknown z_dim_method.');
    end
end

% -------------------------
% Group Z into bins
% -------------------------
group_id = make_groups(Z_group_base, prm.nbins, prm.group_method);
ug = unique(group_id(~isnan(group_id)));
K = numel(ug);

% -------------------------
% Decide plotting mode
% -------------------------
if strcmpi(prm.mode, 'auto')
    if size(Xp,2) == 1
        mode_use = 'univariate';
    else
        mode_use = 'joint';
    end
else
    mode_use = lower(prm.mode);
end

% -------------------------
% Prepare X for plotting
% -------------------------
pcaX = [];
switch mode_use
    case 'univariate'
        if size(Xp,2) == 1
            X_plot = Xp(:,1);
            x_label = 'X';
        else
            switch lower(prm.x_method)
                case 'direct'
                    X_plot = Xp(:, prm.x_dims(1));
                    x_label = sprintf('X_%d', prm.x_dims(1));
                case 'pca'
                    [coeffX, scoreX, latentX, ~, explainedX] = pca(Xp);
                    X_plot = scoreX(:,1);
                    x_label = 'PC1(X)';
                    pcaX.coeff = coeffX;
                    pcaX.score = scoreX;
                    pcaX.latent = latentX;
                    pcaX.explained = explainedX;
                otherwise
                    error('Unknown x_method.');
            end
        end

    case 'joint'
        if size(Xp,2) < 2 && ~strcmpi(prm.x_method,'pca')
            error('Joint mode requires X to have at least 2 dimensions.');
        end

        switch lower(prm.x_method)
            case 'direct'
                if numel(prm.x_dims) ~= 2
                    error('x_dims must have two entries in joint mode.');
                end
                X_plot = Xp(:, prm.x_dims);
                x_label = {sprintf('X_%d', prm.x_dims(1)), sprintf('X_%d', prm.x_dims(2))};

            case 'pca'
                [coeffX, scoreX, latentX, ~, explainedX] = pca(Xp);
                X_plot = scoreX(:,1:2);
                x_label = {'PC1(X)', 'PC2(X)'};
                pcaX.coeff = coeffX;
                pcaX.score = scoreX;
                pcaX.latent = latentX;
                pcaX.explained = explainedX;

            otherwise
                error('Unknown x_method.');
        end

    otherwise
        error('Unknown mode.');
end

% title helper
if isempty(prm.title_str)
    title_main = sprintf('Joint structure | X: %s space | Z: %s space', ...
        x_space_str, z_space_str);
else
    title_main = sprintf('%s | X: %s space | Z: %s space', ...
        prm.title_str, x_space_str, z_space_str);
end

% -------------------------
% Plot
% -------------------------
switch mode_use
    case 'univariate'
        figure;
        switch lower(prm.plot_type_1d)
            case 'hist'
                tiledlayout(K,1, 'Padding','compact', 'TileSpacing','compact');
                for i = 1:K
                    nexttile;
                    idx = group_id == ug(i);
                    histogram(X_plot(idx), 'Normalization', 'pdf');
                    xlabel(x_label);
                    ylabel('PDF');
                    title(sprintf('Group %d | X:%s Z:%s', ug(i), x_space_str, z_space_str));
                end

            case 'box'
                boxchart(categorical(group_id), X_plot);
                xlabel(sprintf('Z group (%s space)', z_space_str));
                ylabel(sprintf('%s (%s space)', x_label, x_space_str));
                title(title_main);

            otherwise
                error('Unknown plot_type_1d.');
        end

    case 'joint'
        switch lower(prm.plot_type_2d)
            case 'overlay_scatter'
                figure; hold on;
                cmap = lines(K);

                h_ellipse = gobjects(K,1);

                for i = 1:K
                    idx = group_id == ug(i);
                    Xi = X_plot(idx, :);

                    scatter(Xi(:,1), Xi(:,2), 18, ...
                        'MarkerFaceColor', cmap(i,:), ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerFaceAlpha', 0.25);

                    if size(Xi,1) > 5
                        h_ellipse(i) = plot_cov_ellipse(Xi, cmap(i,:), 2);
                    end
                end

                h_overall = plot_cov_ellipse(X_plot, [0 0 0], 2, '--', 2);

                xlabel(sprintf('%s (%s)', x_label{1}, x_space_str));
                ylabel(sprintf('%s (%s)', x_label{2}, x_space_str));
                title(title_main);
                axis equal;
                box on;

                lgd_txt = arrayfun(@(x) sprintf('Group %d', x), ug, 'UniformOutput', false);
                lgd_txt{end+1} = 'Overall';
                legend([h_ellipse; h_overall], lgd_txt, 'Location', 'best');

            case 'subplot_scatter'
                figure;
                tiledlayout(1,K, 'Padding','compact', 'TileSpacing','compact');
                for i = 1:K
                    nexttile;
                    idx = group_id == ug(i);
                    scatter(X_plot(idx,1), X_plot(idx,2), 15, 'filled', ...
                        'MarkerFaceAlpha', 0.5);
                    xlabel(sprintf('%s (%s)', x_label{1}, x_space_str));
                    ylabel(sprintf('%s (%s)', x_label{2}, x_space_str));
                    title(sprintf('Group %d | Z:%s', ug(i), z_space_str));
                end

            case 'contour'
                figure;
                tiledlayout(1,K, 'Padding','compact', 'TileSpacing','compact');
                for i = 1:K
                    nexttile;
                    idx = group_id == ug(i);
                    plot_density2d(X_plot(idx,1), X_plot(idx,2));
                    xlabel(sprintf('%s (%s)', x_label{1}, x_space_str));
                    ylabel(sprintf('%s (%s)', x_label{2}, x_space_str));
                    title(sprintf('Group %d', ug(i)));
                end

            otherwise
                error('Unknown plot_type_2d.');
        end
end

% -------------------------
% Output
% -------------------------
out = struct();
out.group_id = group_id;
out.Z_group_base = Z_group_base;
out.X_plot = X_plot;
out.X_used = Xp;
out.Z_used = Zp;
out.X_space = x_space_str;
out.Z_space = z_space_str;
out.params = prm;
out.pcaX = pcaX;
out.pcaZ = pcaZ;
out.valid_idx = valid_idx;

end




















