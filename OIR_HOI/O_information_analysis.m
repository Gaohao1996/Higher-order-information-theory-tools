close all; 
clear;
clc;

%% data import
%Oscillation types 
% - 1: Delta: 0.5-4 Hz (2-4Hz; avoid volume conduction)
% - 2: Theta: 4-8 Hz 
% - 3: Alpha: 8-13 Hz 
% - 4: Beta: 13 - 30 Hz 
% - 5: Gamma: 30 - 80 Hz 
% - 6: HFO: 80 - 200 Hz
period = 'ict';    % 或 'pre'
bands = [2 4; 4 8; 8 13; 13 30; 30 80; 80 200];
nbands = size(bands,1);
OIR_results_human1 = load(sprintf('%s_humandata_1secwin.mat', period)).OIR_results_human;
OIR_results_human2 = load(sprintf('%s_humandata_2secwin.mat', period)).OIR_results_human;
OIR_results_human3 = load(sprintf('%s_humandata_3secwin.mat', period)).OIR_results_human;
OIR_results_human4 = load(sprintf('%s_humandata_4secwin.mat', period)).OIR_results_human;
Oscillation_types = {'Delta','Theta','Alpha','Beta','Gamma','HFO'};
information_flows = {'total information'; 'directed information';'reverse information';'couplings'};


%% O-information
O_spectrum = struct();
O_spectrum(1).total = integrate_OIR_band(OIR_results_human1, []);
O_spectrum(2).total  = integrate_OIR_band(OIR_results_human2, []);
O_spectrum(3).total = integrate_OIR_band(OIR_results_human3, []);
O_spectrum(4).total = integrate_OIR_band(OIR_results_human4, []);
% O_spectrum1 = integrate_OIR_band(OIR_results_human1, []);
% O_spectrum2 = integrate_OIR_band(OIR_results_human2, []);
% O_spectrum3 = integrate_OIR_band(OIR_results_human3, []);
% O_spectrum4 = integrate_OIR_band(OIR_results_human4, []);

%% normalize the O-information by time
for i =1:4
    O_spectrum(i).total.O12= O_spectrum(i).total.O12/i;
end

% colors = [
%     44 160 44    % green
%     31 119 180   % blue
%     214 39 40    % red
%     255 127 14   % orange
%     ]/255;
% 
% figure;
%     for i =1:4
%         ax = subplot(2,2,i);
%         [~, ~] = boxplot_diy(O_spectrum(i).total.O12, ...
%             'Parent', ax, ...
%             'YLabel', sprintf('O-information (average by time)'), ...
%             'Title',  sprintf('O-information accross trials with window length %ds',i), ...
%             'Colors', colors(i,:), ...         % 设置每个窗口的边线色
%             'BoxFaceAlpha', 0.18, ...     % 给箱体加半透明填充
%             'Symbol','k+' );              % 离群点标记（可设 '' 隐藏）
%     end 
%     

for b =1:nbands%nbands
    Oscillation_index = b;
    Band_width = bands(Oscillation_index,2) - bands(Oscillation_index,1);
    %% O-information calculation (integration)
    O_spectrum3 = integrate_OIR_band_new(OIR_results_human3.result, OIR_results_human3.f, bands(Oscillation_index,:), ...
    'BandwidthNormalize', true);
%     (OIR_results_human3, bands(Oscillation_index,:));
    %% Boxplot (one band)
    O_spectrum_norm = zeros(size(O_spectrum3.O12,1),size(O_spectrum3.O12,2),4);
    O_spectrum_norm(:,:,1) = O_spectrum3.O12;
    O_spectrum_norm(:,:,2) = O_spectrum3.O1_2;
    O_spectrum_norm(:,:,3) = O_spectrum3.O2_1;
    O_spectrum_norm(:,:,4) = O_spectrum3.O1o2;
%     /3/Band_width;
    colors = [
            44 160 44    % green
            31 119 180   % blue
            214 39 40    % red
            255 127 14   % orange    
    ]/255;
    % title_subplot = sprint
    figure;
    for i =1:4
        ax = subplot(2,2,i);
        [~, ~] = boxplot_diy(O_spectrum_norm(:,:,i), ...
            'Parent', ax, ...
            'YLabel', sprintf('O-information'), ...
            'Title',  sprintf('O-information(%s) accross trials(%s) over time windows in %s band', period, information_flows{i},Oscillation_types{Oscillation_index}), ...
            'Colors', colors(i,:), ...         % 设置每个窗口的边线色
            'BoxFaceAlpha', 0.18, ...     % 给箱体加半透明填
            'Symbol','k+', 'YScaleFactor', 1e3, ...
                'AutoTight', true);              % 离群点标记（可设 '' 隐藏）
    end 
end


% O-information accross trials by shifting time windows
%% Boxplot with comparision among different time windows(resampled to fixed K windows)
% mats = {O_spectrum(1).total.O12, O_spectrum(2).total.O12 ,O_spectrum(3).total.O12,O_spectrum(4).total.O12};
% labels = {'1s win', '2s win', '3s win', '4s win'};     
% Ks = [4 5 6];
% OI_renomalized_5 = zeros(size(O_spectrum(1).total,1),5,numel(mats));
% % OI_renomalized_4 = zeros(size(O12_spectrum1_norm,1),4,numel(mats));
% 
% for Ki = 1:numel(Ks)
%     K = Ks(Ki);
%     Renormalize_OI = grouped_boxplot_resampled(mats, labels, K, ...
%         'Interp','linear', ...
%         'YLabel','O-information (integral)', ...
%         'Title','O-information with different window length', ...
%         'Test','none', ...
%         'FDR',true, ...
%         'Alpha',0.05);
%     summarize_stability(Renormalize_OI,labels)
%     
% %     if K == 5
% %         for i = 1:numel(mats)
% %             OI_renomalized_5(:,:,i) = Renormalize_OI{i}(:,:);
% %         end     
% %     end
% end


%% significant test and effect
% % K =5
% T3_4 = signrank_bins_paired(OI_renomalized_5(:,:,3), OI_renomalized_5(:,:,4), '3s','4s', 0.05, true);
% T3_1 = signrank_bins_paired(OI_renomalized_5(:,:,3), OI_renomalized_5(:,:,1), '3s','1s', 0.05, true);
% T3_2 = signrank_bins_paired(OI_renomalized_5(:,:,3), OI_renomalized_5(:,:,2), '3s','2s', 0.05, true);





