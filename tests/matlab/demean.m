function x=demean(x,dim)
% Source: https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m
% Original upstream function, unchanged except these source comments.

% DEMEAN(X) 
% Removes the Average or mean value.
%
% DEMEAN(X,DIM)
% Removes the mean along the dimension DIM of X. 

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
