% Independent Octave/MATLAB-compatible checks of the four WEiDA stages.
% Run with working directory containing input.mat. Packages only for Octave.
if exist('OCTAVE_VERSION','builtin')
  pkg load signal; pkg load statistics;
end
load('input.mat');
rng(284, 'twister');
[b,a]=butter(2,[.04 .07]/(1/(2*tr)));
filtered=zeros(size(bold)); phase=zeros(size(bold));
for region=1:size(bold,1)
  ts=detrend(bold(region,:)); ts=ts-mean(ts);
  filtered(region,:)=filtfilt(b,a,ts);
  phase(region,:)=angle(hilbert(filtered(region,:)));
end
native_features=zeros(size(bold,2),size(bold,1));
% Same shared phase isolates eigensolver behavior from filter differences.
shared_features=native_features;
for t=1:size(bold,2)
  for mode=1:2
    if mode==1; phi=phase(:,t); else; phi=shared_phase(:,t); end
    matrix=cos(phi-phi');
    [V,L]=eigs(matrix,6,'lm');
    v=V*diag(L)/sum(diag(L));
    if mean(v>0)>.5 || (mean(v>0)==.5 && sum(v(v>0))>-sum(v(v<0))); v=-v; end
    if mode==1; native_features(t,:)=v'; end
    if mode==2; shared_features(t,:)=v'; end
  end
end
% Shared separated clusters and initial centers isolate K-means semantics.
[idx,C,sumd,D]=kmeans(cluster_input,4,'start',initial_centers,'replicates',1,'maxiter',500,'distance','sqeuclidean');
ss=silhouette(cluster_input,idx,'sqeuclidean');
% Shared prescribed labels isolate stages 3 and 4 from earlier differences.
K=4; n=length(lengths); P=zeros(n,K); LT=nan(n,K); transitions=zeros(n,K,K);
prob=zeros(1,K);
for c=1:K; prob(c)=mean(labels==c); end
[~,order]=sort(prob,'descend');
valid=zeros(n,K); offset=0;
for s=1:n
  seq=labels(offset+1:offset+lengths(s)); offset=offset+lengths(s);
  for c=1:K
    P(s,c)=mean(seq==c);
    active=seq==c; starts=find(diff(active)==1); ends=find(diff(active)==-1);
    if length(ends)>length(starts); ends(1)=[];
    elseif length(starts)>length(ends); starts(end)=[];
    elseif ~isempty(starts) && ~isempty(ends) && starts(1)>ends(1); ends(1)=[]; starts(end)=[];
    end
    if ~isempty(starts) && ~isempty(ends); LT(s,c)=mean(ends-starts)*tr; end
  end
  for i=1:K
    denom=sum(seq(1:end-1)==order(i)); valid(s,i)=denom>0;
    if denom>0
      for j=1:K
        transitions(s,i,j)=sum(seq(1:end-1)==order(i) & seq(2:end)==order(j))/denom;
      end
    end
  end
end
group=zeros(2,K,K);
for g=1:2
  for i=1:K
    eligible=find(conditions==g & valid(:,i));
    for j=1:K
      if isempty(eligible); group(g,i,j)=NaN;
      else; group(g,i,j)=mean(transitions(eligible,i,j)); end
    end
  end
end
runtime_version=version;
save('reference.mat','b','a','filtered','phase','native_features','shared_features','idx','C','sumd','D','ss','P','LT','transitions','group','order','runtime_version','-v7');
fprintf('Reference saved successfully\n');
