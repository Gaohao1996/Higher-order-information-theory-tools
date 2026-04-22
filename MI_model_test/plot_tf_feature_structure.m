function out = plot_tf_feature_structure(X, Z, mi_fun, fig_title, varargin)
% plot_tf_feature_structure
%
% 输入
% ----
% X        : Nsamp x 1 复数向量（单个时频点/单窗口/单通道的复时频系数）
% Z        : Nsamp x 1 连续变量（外部变量、行为量、病理强度、状态指标等）
% mi_fun   : 互信息函数句柄，格式 I = mi_fun(X, Y)
%            要支持:
%              X: Nsamp x d
%              Y: Nsamp x d2
% fig_title: 图标题
%
% 可选参数
% --------
% 'color_mode'    : 'continuous' 或 'group'
% 'n_groups'      : 分组数（用于 group/contours），默认 3
% 'show_contours' : true/false，默认 true
%
% 输出
% ----
% out 结构体，包含:
%   out.I_spec
%   out.I_amp
%   out.I_phase
%   out.Xri
%   out.Amp
%   out.Ph2D
%   out.Cph2D

p = inputParser;
addParameter(p, 'color_mode', 'continuous');
addParameter(p, 'n_groups', 3);
addParameter(p, 'show_contours', true);
parse(p, varargin{:});

color_mode = p.Results.color_mode;
n_groups   = p.Results.n_groups;
show_contours = p.Results.show_contours;

X = X(:);
Z = Z(:);

if length(X) ~= length(Z)
    error('X 和 Z 长度必须一致');
end

Nsamp = length(X);

%% --- three representations ---
Amp = abs(X);
Ph  = angle(X);

Xri   = [real(X), imag(X)];
Ph2D  = [cos(Ph), sin(Ph)];

% normalized spectrum on unit circle
Xphs = X ./ abs(X);

% copula transformed phase representation
Cph2D = [copnorm(real(Xphs)), copnorm(imag(Xphs))];

%% --- MI computation ---
I_spec  = mi_fun(Xri, Z);
I_amp   = mi_fun(Amp, Z);
I_phase = mi_fun(Ph2D, Z);

%% --- grouping for visualization ---
%other optional method:equal;
[group_id, group_edges] = make_groups_from_Z(Z, n_groups);

%% --- figure ---
figure('Color', 'w', 'Position', [100 80 1500 700]);

% =========================================================
% 1. Original complex spectrum
% =========================================================
subplot(2,3,1); hold on;
plot_complex_scatter(real(X), imag(X), Z, group_id, color_mode);
xlabel('Real');
ylabel('Imag');
title('Original complex spectrum');
axis square;
grid on;
box off;

% =========================================================
% 2. Amplitude vs Z
% =========================================================
subplot(2,3,2); hold on;
switch lower(color_mode)
    case 'continuous'
        scatter(Z, Amp, 14, Z, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
        colormap(parula);
        colorbar;
    case 'group'
        cols = lines(n_groups);
        for g = 1:n_groups
            idx = group_id == g;
            scatter(Z(idx), Amp(idx), 14, cols(g,:), 'filled', ...
                'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
        end
end
xlabel('Z');
ylabel('Amplitude');
title('Amplitude');
grid on;
box off;

% =========================================================
% 3. Phase / direction on unit circle
% =========================================================
subplot(2,3,3); hold on;
plot_complex_scatter(real(Xphs), imag(Xphs), Z, group_id, color_mode);
th = linspace(0, 2*pi, 400);
plot(cos(th), sin(th), 'k-', 'LineWidth', 1);
xlabel('Real');
ylabel('Imag');
title('Phase / direction (normalized spectrum)');
axis square;
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
grid on;
box off;

% =========================================================
% 4. Copula transformed complex spectrum (optional but useful)
% =========================================================
subplot(2,3,4); hold on;
CXri = [copnorm(real(X)), copnorm(imag(X))];
plot_complex_scatter(CXri(:,1), CXri(:,2), Z, group_id, color_mode);

if show_contours
    draw_group_contours(CXri, group_id, n_groups);
end

xlabel('copnorm(Re)');
ylabel('copnorm(Im)');
title('Copula-transformed spectrum');
axis square;
grid on;
box off;

% =========================================================
% 5. Copula transformed phase
% =========================================================
subplot(2,3,5); hold on;
plot_complex_scatter(Cph2D(:,1), Cph2D(:,2), Z, group_id, color_mode);

if show_contours
    draw_group_contours(Cph2D, group_id, n_groups);
end

xlabel('copnorm(cos\phi)');
ylabel('copnorm(sin\phi)');
title('Copula-transformed phase');
axis square;
grid on;
box off;

% =========================================================
% 6. MI bar plot
% =========================================================
subplot(2,3,6);
bar([I_spec, I_amp, I_phase], 0.65);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Spectrum','Amplitude','Phase'});
ylabel('MI (bits)');
title('MI comparison');
box off;
grid on;

sgtitle(fig_title, 'FontSize', 16, 'FontWeight', 'bold');

%% --- outputs ---
out = struct();
out.I_spec = I_spec;
out.I_amp = I_amp;
out.I_phase = I_phase;
out.Xri = Xri;
out.Amp = Amp;
out.Ph2D = Ph2D;
out.Cph2D = Cph2D;
out.group_edges = group_edges;

% fprintf('\n=== %s ===\n', fig_title);
% fprintf('I_spec  = %.6f bits\n', I_spec);
% fprintf('I_amp   = %.6f bits\n', I_amp);
% fprintf('I_phase = %.6f bits\n', I_phase);

end