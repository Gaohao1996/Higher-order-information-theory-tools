function [gc_red_surr, gc_syn_surr] = CTE_surr(x, p, d, path_red, path_syn, nsurr, source_index, target_index,delayZ, phase)
% CTE_surr
% surrogate test for redundancy/synergy conditional paths
%
% Inputs:
%   x            : T x M data
%   p, d         : parameters for CTE_iteration
%   path_red     : condition-variable indices only, e.g. [c1 c2 ...]
%   path_syn     : condition-variable indices only, e.g. [c1 c2 ...]
%   nsurr        : number of surrogates
%   source_index : source variable index in original x
%   target_index : target variable index in original x
%   phase        : 1 -> circular shift; otherwise IAAFT
%
% Outputs:
%   gc_red_surr  : nsurr x length(path_red)
%   gc_syn_surr  : nsurr x length(path_syn)
%
% Notes:
%   For the dr-th condition variable, the tested model is:
%       [source, target, path(1), ..., path(dr-1), surrogate(path(dr))]
%   In this reduced matrix, source=1 and target=2.

    nRed = length(path_red);
    nSyn = length(path_syn);

    gc_red_surr = NaN(nsurr, nRed);
    gc_syn_surr = NaN(nsurr, nSyn);

    % ---------- redundancy path ----------
    for dr = 1:nRed
        prev_cond = path_red(1:dr-1);

        if isempty(prev_cond)
            xx = x(:, [source_index, target_index]);
        else
            xx = x(:, [source_index, target_index, prev_cond]);
        end

        for k = 1:nsurr
            % surrogate current tested condition variable
            if phase == 1
                N = size(x,1);
                shift = randi([round(0.1*N), round(0.9*N)]);
                ys = circshift(x(:, path_red(dr)), shift);
            else
                ys = surr_iaafft(x(:, path_red(dr)));
            end

            xxx = [xx, ys];

            % in xxx: col1=source, col2=target
            gc = CTE_iteration(xxx, p, d, 1, 2, delayZ, phase);
            gc_red_surr(k, dr) = mean(gc);
        end
    end

    % ---------- synergy path ----------
    for dr = 1:nSyn
        prev_cond = path_syn(1:dr-1);

        if isempty(prev_cond)
            xx = x(:, [source_index, target_index]);
        else
            xx = x(:, [source_index, target_index, prev_cond]);
        end

        for k = 1:nsurr
            if phase == 1
                N = size(x,1);
                shift = randi([round(0.1*N), round(0.9*N)]);
                ys = circshift(x(:, path_syn(dr)), shift);
            else
                ys = surr_iaafft(x(:, path_syn(dr)));
            end

            xxx = [xx, ys];

            % in xxx: col1=source, col2=target
            gc = CTE_iteration(xxx, p, d, 1, 2,delayZ, phase);
            gc_syn_surr(k, dr) = mean(gc);
        end
    end
end



% function [gc_red_surr,gc_syn_surr]=CTE_surr(x,p,d, drivers_red,drivers_syn,nsurr, souce_index, target_index,condition_index,phase) 
% %i -> j
% gc_red_surr = zeros(nsurr,size(x,2));
% gc_syn_surr = zeros(nsurr,size(x,2));
% for dr=condition_index% 3:kmax+2
%     xx=x(:,drivers_red(1:dr-1));
%     for k=1:nsurr
%         if phase == 1
%             N = size(x,1);
%             shift = randi([round(0.1*N), round(0.9*N)]);
%             ys = circshift(x(:,drivers_red(dr)), shift);
%         else
%             ys = surr_iaafft(x(:,drivers_red(dr)));
%         end
% %         disp(std(ys))  % 看一下是否是常数
%         xx(:,dr)=ys;
%         xxx=xx;
%         %1）Gaussian Copula 
%         gc=CTE_iteration(xxx,p,d,1,2,phase);
%         %2）Linear Causal filter 
% %          gc=LCF_TE_multivariate(xxx,p,d,souce_index,target_index,[]);
% %         disp(gc);
%         Gs=mean(gc);
%         gc_red_surr(k,dr)=Gs;
% %         disp(gc_red_surr(k,dr));
%     end
% end
% 
% for dr=condition_index %3:kmax+2
%     xx=x(:,drivers_syn(1:dr-1));
%     for k=1:nsurr
%          if phase == 1
%             N = size(x,1);
%             shift = randi([round(0.1*N), round(0.9*N)]);
%             ys = circshift(x(:,drivers_syn(dr)), shift);
%         else
%             ys = surr_iaafft(x(:,drivers_syn(dr)));
%          end 
%         xx(:,dr)=ys;
%         xxx=xx;
%         %1） Gaussian Copula       
%         gc=CTE_iteration(xxx,p,d,1,2,phase);
%         %2)   Linear Causal filter 
% %         gc=LCF_TE_multivariate(xxx,p,d,souce_index,target_index,[]);
%         Gs=mean(gc);
%         gc_syn_surr(k,dr)=Gs;
%     end
% end
