function PTRANSITION = original_transitions(IDX, NumClusters)
% 出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/WEiDA_fix_cluster.m
% 来自末尾转移概率段；输入单次扫描的状态标签，原算法语句不变。
PTRANSITION=zeros(NumClusters,NumClusters);
i=1;
for c1=1:NumClusters
    j=1;
    for c2=1:NumClusters
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1)));
        end
        j=j+1;
    end
    i=i+1;
end
end
