Data = load('Matrix_2D.mat'); %(group5_BI_beforetest)
EEG_all= Data.Averaged_Samples;
% EEG_cond = your_data_cond2;   % shape: N x P
fs = 200;                      % sample_rate

figure;
%% power spectrum
for i =1:4
x = EEG_all(:,i);
N = size(x,1);
X = fft(x, N);
P2 = (abs(X).^2) / (fs*N);        % 双边功率谱密度
P1 = P2(1:N/2);                   % 单边
P1(2:end-1) = 2*P1(2:end-1);

f = (0:(N/2-1)) * (fs/N);
subplot(4,1,i)
plot(f, P1);
hold on;
xlabel('Frequency (Hz)');
ylabel('PSD (V^2/Hz)');
title(sprintf('Power Spectrum (Channel %d) at sampling rate %d Hz', i, fs));
end



% win_sec = 10;                   % 1–4 s 皆可
% win = hamming(round(win_sec*fs));
% noverlap = round(0.5*length(win));
% nfft = 2^nextpow2(length(win));
% 
% [Pxx, f] = pwelch(x, win, noverlap, nfft, fs, 'onesided');  % V^2/Hz
% 
% figure;
% plot(f, Pxx); xlim([0 100]);
% xlabel('Frequency (Hz)'); ylabel('PSD (V^2/Hz)'); grid on

