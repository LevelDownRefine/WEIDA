function OPs = original_metastability(Weighted_Eig, state_ts, Vemp)
% 出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/MSanalysis.m
% 来自 DBS 段的状态离散度公式；只把受试者维度和固定 4 状态改为函数参数。
OPs=zeros(1,size(Vemp,2));
for st=1:size(Vemp,2)
    OPs(st)=mean(sum((Weighted_Eig(state_ts==st,:)-Vemp(:,st)').^2) / (sum(state_ts==st)-1));
end
end
