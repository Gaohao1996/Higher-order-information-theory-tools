this_dir = fileparts(mfilename('fullpath'));   % 当前脚本所在文件夹
project_root = fileparts(this_dir);            % 上一级：project_code
addpath(genpath(fullfile(project_root,'HOIs_tools')));

close all;
clear;
clc;

N = 2000;
Z = randn(N,1);   % continuous source variable
model_id = 2;     % 1 phase-only, 2 amp-only, 3 both, 4 mixed
A = zeros(N,1);
phi = zeros(N,1);

sigmaA = 0.3;
sigmaPhi = 0.3;

for i = 1:N
    z = Z(i);

    switch model_id
        case 1   % phase-only
            A(i) = chi2rnd(6);
            phi(i) = mod(1.2*z + sigmaPhi*randn, 2*pi);

        case 2   % amplitude-only
            A(i) = exp(0.8*z + sigmaA*randn);
            phi(i) = -pi + 2*pi*rand; %uniform 
            % phi(i) = circ_vmrnd(0, 3, 1);  %von mises

        case 3   % both amplitude and phase
            A(i) = exp(0.8*z + sigmaA*randn);
            phi(i) = mod(1.2*z + sigmaPhi*randn, 2*pi);

        case 4   % mixed / complementary
            A(i) = exp(0.8*z + sigmaA*randn);
            if z > 0
                phi(i) = mod(0 + sigmaPhi*randn, 2*pi);
            else
                phi(i) = mod(pi + sigmaPhi*randn, 2*pi);
            end
    end
end
X = A .* exp(1i*phi);

Xri  = [real(X), imag(X)];
Amp  = A;
Ph2D = [cos(phi), sin(phi)];

Joint = [Amp, Ph2D];






%% parameter options
% 'color_mode'    : 'continuous' 或 'group'
% 'n_groups'      : 分组数（用于 group/contours），默认 3
% 'show_contours' : true/false，默认 true

n_groups =4; 

titles = {'Phase-only','Amplitude-only','Both'};
out = plot_tf_feature_structure(X, Z, @gcmi_cc, ...
    ['Data structure in complex spectrum domain',
    titles(model_id)],'color_mode', 'continuous','n_groups',n_groups); 


%% KNN/GC mode

% opts = struct('k',6,'zscore',true,'eps',1e-15,'metric','chebychev');
% I_spec  = mi_ksg(Xri,Z );


% I_spec  = gcmi_cc(Xri,Z);
% I_amp   = gcmi_cc(Amp,Z);
% I_phase =  gcmi_cc(Ph2D,Z);
% 
% 
% % bar plot
% figure;
% bar([I_spec, I_amp, I_phase], 0.6);
% set(gca, 'XTick', 1:3, 'XTickLabel', {'Spec(Z;X)','Amp(Z;A)','Phi(Z;Ph2D)'});
% ylabel('MI (bits)');
% title('MI comparison');
% box off;
% ylim([0, max([I_spec, I_amp, I_phase])*1.2 + 0.01]);
% 

