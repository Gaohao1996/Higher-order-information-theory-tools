clc;
clear;



%%Toy_model
% N   : 输出长度（不含 burn-in）
% fs  : 采样率
% f0  : 中心频率，比如 10 Hz
% tau : 耦合延迟（样本点），比如 3~5
% r   : 极点半径，越接近1越窄带，比如 0.98~0.995
N = 6000;
fs =200;
f0 = 10;
tau = 3;
r = 0.985;


x = toy_narrowband_OIR(N, fs, f0, tau, r);

% OIR 专用 opts（避免和你上面 opts 混淆，也更不容易出错）
f_OIR = [];
x_oir = x(:,1:3);   
nfft = 256;
pmax = 10;
fs_oir  = fs;  
opts_oir = struct();
opts_oir.var_order  = 1:3;
opts_oir.Mv         = [1 1 1];
opts_oir.iM         = [1 2 3];       % 关注节点1与2
opts_oir.j          = 3;           % 条件/加入变量节点3
opts_oir.nfft       = nfft;   % 512
opts_oir.pmax       = pmax;   % 10
% opts.bands = [0,5];
opts_oir.cond_names = {'rossler'};

results_oir = oir_flexible({x_oir}, fs_oir, opts_oir);
Roir = results_oir.cond(1);
% 频率轴（只存一次）
if isempty(f_OIR)
    f_OIR = Roir.f;
end



figure;
plot(f_OIR, Roir.dO1_2f, 'LineWidth', 1.6); hold on;             % 方向性：sources → j
plot(f_OIR, Roir.dO1o2f, 'LineWidth', 1.6);                      % 耦合：sources 对 j 的协同/冗余
plot(f_OIR, Roir.dO2_1f, 'LineWidth', 1.6);                      % 反向：j → sources
plot(f_OIR, Roir.dO12f,  '--', 'LineWidth', 1.0);                % （可选）总谱
yline(0,'k:');

xlim([0 20]);
xlabel('Frequency (Hz)'); ylabel('OIR spectral density');
legend('dir (sources→j)','coupling','reverse (j→sources)','total'); 