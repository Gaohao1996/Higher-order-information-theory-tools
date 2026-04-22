close all;
clear;
clc;


index_O = [1,2,3,4];
U_information_total = zeros(7,4); % trials X oscillation types
R_information_total = zeros(7,4);
S_information_total = zeros(7,4);

for Oi_index = index_O
    for trial_index = 1:7
        %% 1.Data import
        file_index = trial_index;
        dataPath ='D:\Matlab_project\dataset\Epilepsy_dataset\seizure8';
        filename1 = sprintf('sz%d_pre_clean.mat', file_index );
        Data = load(fullfile(dataPath, filename1));
        % pre_eeg= load('sz1_pre_clean.mat').pre_eeg_clean;
        eeg_data  = Data.pre_eeg;
        % eeg_alltime = [pre_eeg; ict_eeg]; % [sample (M) X channels (N)]
        
        % Data structure
        fs = 400;  %sample rate
        [M, N] =size(eeg_data);
        Channels = 1:N;
        Samples = 1:M;
        
        %% Parameter setting (time duration, filter, variable index, time window, parmeter optimization, surrogate tests)
        %time duation ( epoch in [1/fs, M/fs])
        start_sample = 1;
        end_sample = M;
        t_start = start_sample/fs;
        t_end = end_sample/fs;
        t_duration = [t_start,t_end];
        
        %% Zero values check in EEG data (necessary for Gaussian copula calculation)
        %Due to finite ADC resolution, the broadband iEEG exhibited discrete amplitude levels
        % (median quantization step ≈ 0.90 in acquisition units), 
        % which led to a high proportion of identical samples after delay embedding. 
        % We applied minimal de-quantization by adding Gaussian noise with σ = 5
        % of the estimated quantization step prior to copula normalization, to ensure valid Gaussian-copula CMI estimation.”
        [eeg_data, rep] = preprocess_ieeg_for_te(eeg_data, 0.05);
        %         oscillation types 
        %oscillation types
        %  - 1: Delta: 2-4 Hz
        %  - 2: Theta: 4-8 Hz
        %  - 3: Alpha: 8-13 Hz
        %  - 4: Beta: 13 - 30 Hz
        %  - 5: Gamma: 30 - 80 Hz
        %  - 6: HFO: 80 - 200 Hz
        Oscillation_index = Oi_index;
        
        % select ROI for TE
        source_channel = [75 76];
        condition_channel = [65 66 67 68 69 70];  %Grid_channels;%setdiff(Grid_channels,EZ_channels);
        target_channel = [1,2,3,4,9,10,11,12,17,18,19,26];
        
        % time window shifting (optional)
        window_length = size(eeg_data,1);   %.g. ceil(0.5*fs)
        overlap_rate = 0;
        
        % Parameters optimization for TE(order and delay)
        % a. order and delay
%         model_orders = 1:5; % not too long for test
%         delays = 50:100;  %based on the required cycles for oscillations
        fc = [3,6,10,20];   
        delays=  floor(0.1*fs/fc(Oscillation_index)):floor(0.6*fs/fc(Oscillation_index));
        if Oscillation_index ==1
                    model_orders = 2:6; % not too long for test
%                     delays = (0.1*400/fc(Oscillation_index)):(0.6*400/fc(Oscillation_index));  %based on the required cycles for oscillations
        elseif Oscillation_index ==2
                    model_orders = 2:5; % not too long for test
%                     delays = 50:100;  %based on the required cycles for oscillations
        elseif Oscillation_index ==3
                    model_orders = 1:4; % not too long for test
%                     delays = 50:100;  %based on the required cycles for oscillations
        else
                    model_orders = 1:3; % not too long for test
%                     delays = 50:100;  %based on the required cycles for oscillations
        end
         
        
        %b. mutual information calculation (optional)
        % Default selection for MI (change the default option in  function named 'gccmi_ccc'  if necessary)
        %biascorrect = true; % whether bias correction should be applied to the esimtated MI
        %demeaned = true;  % already been copula-nomarlized so that no need to change
        %cov = true; % when the covariance matrix is illconditioned use the 'false' button to reduce it (shrinkage matrix)
        
        % TE decomposition and surrogate tests (test amount and p value setting)
        nsurr=99;
        alpha = 0.05;
        
        % Phase mode (1: calculate phase TE ; 0: normal TE)
        phase =1;
        
        
        %% 2. Filter the oscillation (optional, filter parameter need to change in data_filter if sample size changes)
        %default freq band for HF is [80 150] from referrence;
        %select freq band for certain oscillations;
        % Low Frequency band for phase study
        Filtered_EEG = data_filter(eeg_data, fs, t_duration, Channels, Oscillation_index);
        
        
        %% 3.Settting source,target and conditioned variables
        
        %original indexes for source, target and conditioned variables
        source_target_channels = [source_channel,target_channel];
        Channels_idx = [source_channel, target_channel, condition_channel];
        
        %new indexes for source, target and conditioned variables
        sc_num = length(source_channel);
        tc_num = length(target_channel);
        cc_num = length(condition_channel);
        souce_index = 1:sc_num;
        target_index = sc_num+1:sc_num+tc_num;
        condition_index = sc_num+tc_num+1:length(Channels_idx);
        
        %% Extract phase information from the filtered signal
        %phase data from certain band signal
        % Phase_EEG = angle(hilbert(Filtered_EEG)); %Theta Band Signal
        Phase_EEG = angle(hilbert(Filtered_EEG));
        x_all = Phase_EEG(:,Channels_idx);
%         x_all = Filtered_EEG(:,Channels_idx);
        
        %% 4.Time window shifting for dynamical analysis(optional)
        overlap = overlap_rate*window_length; %e.g. 0.5*window_length
        x_windows = time_window_shifting(x_all, window_length, overlap); %channels x samples x time windows
        x_windows = permute(x_windows, [2 1 3]); % samples x channels x time windows
        [~,~,window_number] = size(x_windows);
        % x = x_windows(:,:,1);
        
        
        %% 5.TE decomposition for phase data
        
        %% a). Model parameters optimizing (AIC)
        for window_index = 1: window_number
            x = x_windows(:,:,window_index);
            % AIC for searching best parameter
            [best_p, best_d, ~, ~, ~] = optimize_te_aic_delay_on_Y(x(:,souce_index), x(:,target_index),[], model_orders, delays,phase);
            
            %% b). TE decomposition calculation
            [drivers_red, drivers_syn, g_red, g_syn]=TE_syn_red(x, souce_index, target_index, best_p, best_d,phase);
            
            %% c). Surrogate test
            [g_red_surr,g_syn_surr]=CTE_surr(x ,best_p, best_d, drivers_red,drivers_syn,nsurr,souce_index, target_index,condition_index,phase);
            
            % p value for the randomized test
            p_x = 1:cc_num+1;
            p_values_red = ones(1, cc_num+1);
            p_values_syn = ones(1, cc_num+1);
            for t = 2:cc_num+1
                rank_red = sum(g_red_surr(:, condition_index(t-1)) <= g_red(t));  %
                rank_syn = sum(g_syn_surr(:, condition_index(t-1)) >= g_syn(t));
                p_values_red(t) = (rank_red+1) / (size(g_red_surr, 1) + 1);
                p_values_syn(t) = (rank_syn+1) / (size(g_syn_surr, 1) + 1);
            end
            sig_idx_red = find(p_values_red < alpha);
            sig_idx_syn = find(p_values_syn < alpha);
            
            %% === 基线 ===
            g0 = g_red(1);   % == g_syn(1)
            
            %% === 找“第一次不显著”的步（t），各自独立判定）===
            % 若全显著，则取最后一步
            first_nonsig_red = find(p_values_red(2:end) >= alpha, 1, 'first');
            first_nonsig_syn = find(p_values_syn(2:end) >= alpha, 1, 'first');
            
            if isempty(first_nonsig_red)
                tR_use = numel(g_red);        % 用最后一步
            else
                tR_use = (first_nonsig_red + 1) - 1;   % 用 t-1
            end
            if isempty(first_nonsig_syn)
                tS_use = numel(g_syn);
            else
                tS_use = (first_nonsig_syn + 1) - 1;
            end
            
            % 防御（至少不能小于2；t=2 对应加入第一个变量后的TE）
            tR_use = max(tR_use, 2);
            tS_use = max(tS_use, 2);
            
            gR_use = g_red(tR_use);   % TE(i->j | 冗余前缀)
            gS_use = g_syn(tS_use);   % TE(i->j | 协同前缀)
            
            %% === 信息流分解（简单对比法）===
            R_information = max(0, g0 - gR_use);   % 冗余：基线下降量
            S_information = max(0, gS_use - g0);   % 协同：相对基线的增量
            U_information = max(0, g0 - R_information - S_information); % 剩余独特
  
            U_information_total(trial_index,Oi_index) = U_information;
            R_information_total(trial_index,Oi_index) = R_information;
            S_information_total(trial_index,Oi_index) = S_information;
            
            % 可选：记录使用到的步数便于调试
            rank_red = tR_use;   % 使用到的冗余步（与 g_red 索引一致）
            rank_syn = tS_use;   % 使用到的协同步
            
            
            %     %% 6. plot redundant information change in TE Decomposition
            %     figure
            %     plot(g_red,'-*k');hold on
            %     plot(2:cc_num+1,g_red_surr(:,condition_index),'-or');hold on
            %     plot(p_x (sig_idx_red), g_red(sig_idx_red), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
            %     % xlim([0 3])
            %     for i = 1:length(sig_idx_red)
            %         x_val_red = p_x(sig_idx_red(i));
            %         y_val_red = g_red(sig_idx_red(i));
            %         text(x_val_red, y_val_red, sprintf('%d', Channels_idx(drivers_red(condition_index(sig_idx_red(i)-1)))), ...
            %             'VerticalAlignment', 'top', ...
            %             'HorizontalAlignment', 'center', ...
            %             'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
            %     end
            %     xlabel('Number of Conditioned Variables')
            %     ylabel(sprintf('TE values with order %d and lag %d', best_p, best_d))
            %     title(['Redundant information flow test at time window', num2str(window_index)])
            %     legend('Raw signal','Surrogate Signal')
            %     %% plot synergistic information change in TE Decomposition
            %     figure
            %     plot(g_syn,'-*k');hold on
            %     plot(2:cc_num+1,g_syn_surr(:,condition_index),'-or');
            %     plot(p_x (sig_idx_syn), g_syn(sig_idx_syn), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
            %     % xlim([0 3])
            %     for i = 1:length(sig_idx_syn)
            %         x_val_syn = p_x(sig_idx_syn(i));
            %         y_val_syn = g_syn(sig_idx_syn(i));
            %         text(x_val_syn , y_val_syn, sprintf('%d', Channels_idx(drivers_syn(condition_index(sig_idx_syn(i)-1)))), ...
            %             'VerticalAlignment', 'top', ...
            %             'HorizontalAlignment', 'center', ...
            %             'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
            %     end
            %     xlabel('Number of Conditioned Variables')
            %     ylabel(sprintf('TE values with order %d and lag %d', best_p, best_d))
            %     title(['Synergistic information flow test at time window', num2str(window_index)])
            %     legend('Raw signal','Surrogate Signal')
        end
    end
    
end

