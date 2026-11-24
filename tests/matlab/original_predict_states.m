function IDX = original_predict_states(Weighted_Eig, Center)
% 出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/WEiDA_fix_cluster.m
% 来自第 2 节的最近中心分配；输入固定特征，跳过 BOLD、相位和特征分解。
t_all=size(Weighted_Eig,1); NumClusters=size(Center,1);
IDX=zeros(t_all,1);

for t=1:t_all
    for j=1:NumClusters
        di(j)=sqrt(sum((Weighted_Eig(t,:)-Center(j,:)).^2));  % 固定特征与状态中心的欧氏距离
    end
    [aux indmin]=min(di);   % 按最小欧氏距离确定状态；距离相同则取第一个中心
    IDX(t)=indmin;
end

end
