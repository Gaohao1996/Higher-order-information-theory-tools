dsclose all; 
clear;
clc;

OIR_results_human = struct();
index = 1:7;

%% mode selection: preictal/ictal period (0/1)
mode = 1;
if mode == 0
    period = 'pre';
else
    period = 'ict';
end

%% parameter initialization
fs = 400;                      % sample_rate


%% information path
X1 = [75 76];
X2 = [65 66 67 68 69 70];
X3 =  [1,2,3,4,9,10,11,12,17,18,19,26];
X_channels = [X1, X2, X3];



%% import human epilepsy dataset
for file_index = index
    %% file load
    %% 1) human epilepsy data
    dataPath ='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';
    filename1 = sprintf('sz%d_%s_clean.mat', file_index, period);
    title_name1 = sprintf('OIR spectral density in %s%d', period,file_index );
    Data = load(fullfile(dataPath, filename1));
    
    %% 2) rat epilepsy data
    %dataPath = 'D:\Matlab_project\dataset\Epilepsy_dataset\rat_model\data_used\391_group3_BI';
    % filename2 = sprintf('sz%d_ict_clean.mat', file_index );
    % title_name2 = sprintf('ict%d', file_index );
    % Data = load(fullfile(dataPath, filename2));
   
    
    %% import data matrix(option : multiple task data)
    EEG_all= Data.ict_eeg;
    
    % EEG_all= Data.ict_eeg;
    % EEG_all= Data.Matrix_2D;
    % EEG_all = load('Matrix_2D').Averaged_Samples;
    % EEG_cond = your_data_cond2;   % shape: N x P
       
    % 例子：把通道按 ROI 分组：左顶叶 3通道、右顶叶 4通道、前额 2通道
    % Test X1 + X2 → X3; X3 → X1 + X2
    EEG_all_time = EEG_all(:, X_channels);
    
    %% rat model variables
    % X1 = 2;
    % X2 = [3 4];
    % X3 =  1;
    % X_channels = [X1, X2, X3];
    % EEG_all_time = EEG_all(:, X_channels);
    
    %%  time winodw shifting
    win_sec = 4;
    window_length = ceil(win_sec*fs);%size(EEG_all_time,1);   %e.g. ceil(0.5*fs)
    overlap_rate = 0.5;
    overlap = overlap_rate*window_length;  % e.g. 0.5*window_length
    EEG_windows = time_window_shifting(EEG_all_time, window_length, overlap);  %channels x samples x time windows
    EEG_windows = permute(EEG_windows, [2 1 3]);  % samples x channels x time windows
    [~,~,window_number] = size(EEG_windows);
    
    
    
    %% Analysis in shifting window
    var_order = 1:numel(X_channels); %
    Mv = [numel(X1), numel(X2), numel(X3)]; %(source,mediator,target)
    iM = [1 2 3]; %index of the components
    j  = 3;  % X3 as the added variable
    opts = struct();
    opts.var_order  = var_order;
    opts.Mv         = Mv;
    opts.iM         = iM;
    opts.j          = j;
    opts.nfft       = 512;
    opts.pmax       = 10;          % VAR 选阶上限
    % opts.bands      = {'delta',[1 4]; 'theta',[4 8]; 'alpha',[8 12]; 'beta',[13 30]; 'gamma',[30 80]};
    opts.plot       = false;
    % opts.cond_names = {title_name1 ,'epilepsy seizure'}; %names of the conditions
    opts.cond_names = {title_name1 ,'epilepsy study'}; %names of the conditions
    dO1_2f_all = zeros(opts.nfft, window_number);
    dO1o2f_all = zeros(opts.nfft, window_number);
    dO2_1f_all =  zeros(opts.nfft, window_number);
    dO12f_all = zeros(opts.nfft, window_number);
    
    for i = 1:window_number
        EEG = EEG_windows(:,:,i);
        results = oir_flexible({EEG}, fs, opts);
        %test the error between time and frequency domain
        diag = oir_consistency_diagnostics(results);
%         disp(diag(1).('dO12').relerr)
%          disp(diag(1).('dO1_2').relerr)
%         disp(diag(1).('dO1o2').relerr)
%         disp(diag(1).('dO2_1').relerr)
%         
        if diag(1).('dO1_2').relerr >0.05 || diag(1).('dO1o2').relerr>0.05 || diag(1).('dO2_1').relerr>0.1
            disp('significant error between calculation and theoretical result')
        end
        
        % 取某个带内结果，例如 alpha 带的各分量（条件1）
        R = results.cond(1);
        dO1_2f_all(:,i)  = R.dO1_2f;
        dO1o2f_all(:,i) = R.dO1o2f;
        dO2_1f_all(:,i) = R.dO2_1f;
        dO12f_all(:,i) = R.dO12f;
        
        f = R.f;
        
        % 自动 y 轴范围（含总谱，便于留出上下边距）
        Y = [R.dO1_2f(:); R.dO1o2f(:); R.dO2_1f(:); R.dO12f(:)];
        m = min(Y); M = max(Y); pad = 0.06*(M-m + eps);
        
        
%         figure('Color','w','Name',['OIR spectra - ' results.info.cond_names{1}]);
%         plot(f, R.dO1_2f, 'LineWidth', 1.6); hold on;             % 方向性：sources → j
%         plot(f, R.dO1o2f, 'LineWidth', 1.6);                      % 耦合：sources 对 j 的协同/冗余
%         plot(f, R.dO2_1f, 'LineWidth', 1.6);                      % 反向：j → sources
%         plot(f, R.dO12f,  '--', 'LineWidth', 1.0);                % （可选）总谱
%         yline(0,'k:');
%         
%         xlim([0 200]);
%         ylim([m-pad, M+pad]); grid on;
%         xlabel('Frequency (Hz)'); ylabel('OIR spectral density');
%         legend({'dir (sources→j)','coupling','reverse (j→sources)','total'}, 'Location','best'); legend boxoff;
%         title(sprintf('%s at time window %d', results.info.cond_names{1}, i));
        
        OIR_results_human.result(file_index).win(i).O1_2 = R.dO1_2f;
        OIR_results_human.result(file_index).win(i).O1o2 =R.dO1o2f;
        OIR_results_human.result(file_index).win(i).O2_1 =R.dO2_1f;
        OIR_results_human.result(file_index).win(i).O12 = R.dO12f;
    end
end
OIR_results_human.f = f;
save(sprintf('%s_humandata_%dsecwin.mat', period, win_sec), 'OIR_results_human')




 


