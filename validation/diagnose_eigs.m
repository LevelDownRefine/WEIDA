function diagnose_eigs(folder)
% Black-box diagnostics only. Export numerical inputs/outputs, not MATLAB code.
ref = load(fullfile(folder,'reference.mat'));
indices = [1 2 17 39 77 105 145 190 210];
phases = ref.phase(:,indices);
rng(125,'twister');
phases = [phases, 2*pi*rand(size(phases,1),12)];
% Held-out cases and low-rank/degenerate boundaries; never use them to fit signs.
rng(20261123,'twister');
phases = [phases, 2*pi*rand(size(phases,1),100)];
n = size(phases,1);
boundary = [zeros(n,1), ones(n,1), linspace(0,1e-8,n)', ...
    [zeros(n/2,1); (pi/2)*ones(n/2,1)], (0:n-1)'*2*pi/n];
phases = [phases, boundary];
case_groups = [repmat({'reference'},1,9), repmat({'random_seed_125'},1,12), ...
    repmat({'heldout_seed_20261123'},1,100), ...
    {'constant_zero','constant_one','almost_constant','equal_eigenvalues_blocks', ...
     'equal_eigenvalues_circle'}];
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
    [V,D,flags(i,5)] = eigs(A,6,'largestabs','StartVector',trace(:,2));
    start_values(i,:) = diag(D)'; start_vectors(i,:,:) = reshape(V,[1 n 6]);
end
% Identify the private starting-vector generator by exported numerical behavior.
types = {'mt19937ar','dsfmt19937','mcg16807','mrg32k3a','mlfg6331_64', ...
    'swb2712','shr3cong','threefry4x64_20','philox4x32_10'};
normal_methods = {'Ziggurat','Polar','Inversion'};
stream_trials = struct('type',{},'seed',{},'normal',{},'error',{});
target_start = operator_inputs{1}(:,2);
for j=1:numel(types)
  for seed=[0 1 5489]
    for k=1:numel(normal_methods)
      try
        stream=RandStream(types{j},'Seed',seed,'NormalTransform',normal_methods{k});
        z=randn(stream,n,1); z=z/norm(z);
        trial.type=types{j}; trial.seed=seed; trial.normal=normal_methods{k};
        trial.error=max(abs(z-target_start));
        stream_trials(end+1)=trial;
      catch
        % Unsupported generator/normal-transform combinations are not candidates.
      end
    end
  end
end
% Names only: useful for identifying the actual solver family without copying source.
profile clear; profile on;
eigs(squeeze(matrices(1,:,:)));
info = profile('info'); profile off;
profile_names = {info.FunctionTable.FunctionName};
runtime_version = version;
% Public numerical RNG observations for an independent initialization port.
stream=RandStream('dsfmt19937','Seed',0,'NormalTransform','Ziggurat');
rng_fixture.initial_state=stream.State;
rng_fixture.uniform=rand(stream,2048,1);
stream=RandStream('dsfmt19937','Seed',0,'NormalTransform','Ziggurat');
rng_fixture.normal=randn(stream,2048,1);
rng_fixture.normal_state=stream.State;
save(fullfile(folder,'eigs_diagnosis.mat'),'matrices','native_values','native_vectors', ...
    'repeat_values','repeat_vectors','seed_values','seed_vectors','handle_values', ...
    'handle_vectors','start_values','start_vectors','operator_inputs','flags', ...
    'profile_names','runtime_version','stream_trials','case_groups','rng_fixture','-v7');
end

function [V,D,inputs,flag] = tracedEigs(A)
inputs = zeros(size(A,1),0);
[V,D,flag] = eigs(@multiply,size(A,1),6,'largestabs','IsFunctionSymmetric',true);
    function y = multiply(x)
        inputs = [inputs,x];
        y = A*x;
    end
end
