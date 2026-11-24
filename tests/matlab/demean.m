function x=demean(x,dim)
% 出处： https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m
% 原版上游函数；仅添加出处说明并将注释翻译为中文，算法不变。

% DEMEAN(X) 
% 去除均值。
%
% DEMEAN(X,DIM)
% 沿 X 的 DIM 维度去除均值。 

if(nargin==1),
   dim = 1;
   if(size(x,1) > 1)
      dim = 1;
   elseif(size(x,2) > 1)
      dim = 2;
   end;
end;

dims = size(x);
dimsize = size(x,dim);
dimrep = ones(1,length(dims));
dimrep(dim) = dimsize;

x = x - repmat(mean(x,dim),dimrep);
