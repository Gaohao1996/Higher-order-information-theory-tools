this_dir = fileparts(mfilename('fullpath'));   % 当前脚本所在文件夹
project_root = fileparts(this_dir);            % 上一级：project_code
addpath(genpath(fullfile(project_root,'HOIs_tools')));

% close all;
clear;
clc;

nRep = 20;
allMI = zeros(3,3,nRep);
N = 2000;
sigmaA = 0.3;
sigmaPhi = 0.3;



for r = 1:nRep
    Z = randn(N,1);
    
    for model_id = 1:3
        A = zeros(N,1);
        phi = zeros(N,1);

        for i = 1:N
            z = Z(i);

            switch model_id
                case 1 % phase-only
                    A(i) = chi2rnd(6);
                    phi(i) = mod(1.2*z + sigmaPhi*randn, 2*pi);

                case 2 % amplitude-only
                    A(i) = exp(0.8*z + sigmaA*randn);
                    % phi(i) = -pi + 2*pi*rand;
                    phi(i) = circ_vmrnd(0, 3, 1);  %von mises
                case 3 % mixed 
                    A(i) = exp(0.8*z + sigmaA*randn);
                    phi(i) = mod(1.2*z + sigmaPhi*randn, 2*pi);
            end
        end

        X = A .* exp(1i*phi);
        Xri  = [real(X), imag(X)];
        Amp  = A;    
        phase = mkVar('phi', phi, 'phase');
        % Ph2D = [cos(phi), sin(phi)];
        

        
        %% 1: spectrum; 2:Amplitude; 3:Phase
        
%         opts = struct('k',6,'zscore',true,'eps',1e-15,'metric','chebychev');
%          allMI(model_id,1,r)  = mi_ksg(Xri,Z );
         
        allMI(model_id,1,r) = gcmi_cc(Xri, Z);
        allMI(model_id,2,r) = gcmi_cc(Amp, Z);
        allMI(model_id,3,r) = gcmi_cc(phase, Z);
    end
end

miMean = mean(allMI,3);
miStd  = std(allMI,0,3);

figure('Color','w','Position',[200 200 900 300]);
titles = {'Phase-only','Amplitude-only','Both'};

for k = 1:3
    subplot(1,3,k);
    b = bar(miMean(k,:), 0.6); hold on;
    errorbar(1:3, miMean(k,:), miStd(k,:), 'k.', 'LineWidth', 1.2);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Spec(Z;X)','Amp(Z;A)','Phi(Z;Ph2D)'});
    ylabel('MI (bits)');
    title(titles{k});
    ylim([0, max(miMean(:)+miStd(:))*1.2 + 0.01]);
    box off;
end