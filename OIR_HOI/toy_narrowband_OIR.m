function x = toy_narrowband_OIR(N, fs, f0, tau, r)
% N   : 输出长度（不含 burn-in）
% fs  : 采样率
% f0  : 中心频率，比如 10 Hz
% tau : 耦合延迟（样本点），比如 3~5
% r   : 极点半径，越接近1越窄带，比如 0.98~0.995

if nargin < 1, N = 5000; end
if nargin < 2, fs = 200; end
if nargin < 3, f0 = 10; end
if nargin < 4, tau = 3; end
if nargin < 5, r = 0.985; end

burn = 1000;
T = N + burn + max(tau,2);

% 共振核参数
c = 2 * r * cos(2*pi*f0/fs);
d = -r^2;

% 噪声强度
s1 = 1.0;
s2 = 1.0;
s3 = 0.8;
s4 = 0.8;

% 耦合强度
a13 = 0.35;
a23 = 0.35;
a34 = 0.50;

x = zeros(T,4);
e = randn(T,4);

for t = 3:T
    % X1, X2: 独立窄带源
    x(t,1) = c*x(t-1,1) + d*x(t-2,1) + s1*e(t,1);
    x(t,2) = c*x(t-1,2) + d*x(t-2,2) + s2*e(t,2);

    % X3: 受 X1, X2 共同驱动
    if t-tau >= 1
        drive3 = a13*x(t-tau,1) + a23*x(t-tau,2);
    else
        drive3 = 0;
    end
    x(t,3) = c*x(t-1,3) + d*x(t-2,3) + drive3 + s3*e(t,3);

    % X4: 受 X3 驱动
    if t-tau >= 1
        drive4 = a34*x(t-tau,3);
    else
        drive4 = 0;
    end
    x(t,4) = c*x(t-1,4) + d*x(t-2,4) + drive4 + s4*e(t,4);
end

% 去掉 burn-in
x = x(burn+1:end,:);
end