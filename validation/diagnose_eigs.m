function diagnose_eigs(folder)
% Black-box diagnostics only. Export numerical inputs/outputs, not MATLAB code.
ref = load(fullfile(folder,'reference.mat'));
indices = [1 2 17 39 77 105 145 190 210];
phases = ref.phase(:,indices);
rng(125,'twister');
phases = [phases, 2*pi*rand(size(phases,1),12)];
n = size(phases,1); count = size(phases,2);
matrices = zeros(count,n,n);
native_values = zeros(count,6); native_vectors = zeros(count,n,6);
repeat_values = native_values; repeat_vectors = native_vectors;
seed_values = native_values; seed_vectors = native_vectors;
handle_values = native_values; handle_vectors = native_vectors;
start_values = native_values; start_vectors = native_vectors;
operator_inputs = cell(count,1); flags = zeros(count,5);
for i=1:count
    phi = phases(:,i); A = cos(phi-phi');
    matrices(i,:,:) = reshape(A,[1 n n]);
    [V,D,flags(i,1)] = eigs(A);
    native_values(i,:) = diag(D)'; native_vectors(i,:,:) = reshape(V,[1 n 6]);
    [V,D,flags(i,2)] = eigs(A);
    repeat_values(i,:) = diag(D)'; repeat_vectors(i,:,:) = reshape(V,[1 n 6]);
    rng(104+i,'twister');
    [V,D,flags(i,3)] = eigs(A);
    seed_values(i,:) = diag(D)'; seed_vectors(i,:,:) = reshape(V,[1 n 6]);
    [V,D,trace,flags(i,4)] = tracedEigs(A);
    handle_values(i,:) = diag(D)'; handle_vectors(i,:,:) = reshape(V,[1 n 6]);
    operator_inputs{i} = trace;
    [V,D,flags(i,5)] = eigs(A,6,'largestabs','StartVector',trace(:,1));
    start_values(i,:) = diag(D)'; start_vectors(i,:,:) = reshape(V,[1 n 6]);
end
% Names only: useful for identifying the actual solver family without copying source.
profile clear; profile on;
eigs(squeeze(matrices(1,:,:)));
info = profile('info'); profile off;
profile_names = {info.FunctionTable.FunctionName};
runtime_version = version;
save(fullfile(folder,'eigs_diagnosis.mat'),'matrices','native_values','native_vectors', ...
    'repeat_values','repeat_vectors','seed_values','seed_vectors','handle_values', ...
    'handle_vectors','start_values','start_vectors','operator_inputs','flags', ...
    'profile_names','runtime_version','-v7');
end

function [V,D,inputs,flag] = tracedEigs(A)
inputs = zeros(size(A,1),0);
[V,D,flag] = eigs(@multiply,size(A,1),6,'largestabs','IsFunctionSymmetric',true);
    function y = multiply(x)
        inputs = [inputs,x];
        y = A*x;
    end
end
