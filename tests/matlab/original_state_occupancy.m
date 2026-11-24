function [P, LT] = original_state_occupancy(Ctime, K, TR)
% 出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m
% 来自第 3 节；仅用函数参数替代外层受试者/条件索引，保留首尾截断段规则。
P=zeros(1,K); LT=nan(1,K);
for c=1:K
                P(c)=mean(Ctime==c);

                % 平均持续时间
                Ctime_bin=Ctime==c;

                % 查找进入与离开该状态的位置
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);

                % 排除从扫描起点开始或延伸至扫描终点的不完整状态段
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=NaN;
                end
                LT(c)=mean(C_Durations)*TR;
end
end
