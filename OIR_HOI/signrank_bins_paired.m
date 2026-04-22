function T = signrank_bins_paired(ReA, ReB, labelA, labelB, alpha, doFDR)
% ReA/ReB: trials×K（同一 trials 次序），按秒归一化且重采样到同一 K
if nargin<5, alpha=0.05; end
if nargin<6, doFDR=true; end

K  = size(ReA,2);
Bin = (1:K).'; n = nan(K,1); p = n; pFDR = n; sig = false(K,1);
medA=n; medB=n; iqrA=n; iqrB=n; HL=n;

for k=1:K
    a = ReA(:,k); b = ReB(:,k);
    msk = ~(isnan(a)|isnan(b));
    a = a(msk); b = b(msk);
    n(k) = numel(a);
    if n(k)<3, p(k)=NaN; medA(k)=NaN; medB(k)=NaN; iqrA(k)=NaN; iqrB(k)=NaN; HL(k)=NaN; continue; end

    medA(k)=median(a); medB(k)=median(b);
    iqrA(k)=diff(quantile(a,[.25 .75])); 
    iqrB(k)=diff(quantile(b,[.25 .75]));
    p(k) = signrank(a,b,'method','exact');       % 配对 exact
    HL(k)= median(a-b);                          % Hodges–Lehmann 配对中位数差
end

pFDR = p;
if doFDR
    pFDR = bh_fdr_local(p);
end
sig = pFDR<alpha;

T = table(Bin,n,medA,iqrA,medB,iqrB,HL,p,pFDR,sig, ...
    'VariableNames',{'Bin','n',[labelA '_med'],[labelA '_IQR'],[labelB '_med'],[labelB '_IQR'], ...
                     'HL_AminusB','p','p_FDR','sig_FDR'});
end

function p_adj = bh_fdr_local(p)
    p_adj = nan(size(p)); v=~isnan(p); if ~any(v), return; end
    [ps,idx]=sort(p(v)); m=numel(ps); q=ps.*m./(1:m);
    for i=m-1:-1:1, q(i)=min(q(i),q(i+1)); end
    tmp=nan(size(p)); tmp(v)=q(invperm_local(idx)); p_adj=tmp;
end
function y=invperm_local(idx); y=zeros(size(idx)); y(idx)=1:numel(idx); end

