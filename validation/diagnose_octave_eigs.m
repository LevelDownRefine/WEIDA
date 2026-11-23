function diagnose_octave_eigs(folder)
% Compare engines using literal MATLAB matrices, without regenerating phases.
if exist('OCTAVE_VERSION','builtin') == 0
    error('Run this diagnostic in Octave.');
end
ref = load(fullfile(folder,'eigs_diagnosis.mat'));
rng(284,'twister');
features = zeros(size(ref.matrices,1),size(ref.matrices,2));
flags = zeros(size(ref.matrices,1),1);
for i=1:size(ref.matrices,1)
    [V,D,flags(i)] = eigs(squeeze(ref.matrices(i,:,:)));
    v = V*diag(D)/sum(diag(D));
    if mean(v>0)>.5 || (mean(v>0)==.5 && sum(v(v>0))>-sum(v(v<0)))
        v = -v;
    end
    features(i,:) = v';
end
runtime_version = OCTAVE_VERSION;
save(fullfile(folder,'octave_eigs_diagnosis.mat'),'features','flags','runtime_version','-v7');
end
