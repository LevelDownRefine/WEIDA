function export_weida_features(folder)
% Shared Python preprocessing has already run. Preserve native eigs signs.
if nargin < 1; folder = fileparts(mfilename('fullpath')); end
if exist('OCTAVE_VERSION','builtin')
    pkg load signal;
    runtime_engine = 'Octave';
else
    runtime_engine = 'MATLAB';
end
output = fullfile(folder,'features.mat');
assert(~isfile(output),'features.mat already exists');
input = load(fullfile(folder,'input.mat'));
lengths = cellfun(@(x) size(x,2),input.scans);
n_regions = numel(input.regions);
vectors = zeros(sum(lengths),n_regions);
offset = 0;
for s=1:numel(input.scans)
    amplitudes = input.scans{s};
    assert(size(amplitudes,1)==n_regions && all(isfinite(amplitudes(:))));
    phase = zeros(size(amplitudes));
    for region=1:n_regions
        phase(region,:) = angle(hilbert(amplitudes(region,:)));
    end
    for t=1:lengths(s)
        phi = phase(:,t);
        [V,D] = eigs(cos(phi-phi'));
        v = zeros(n_regions,1);
        for i=1:size(V,2); v = v+D(i,i)*V(:,i); end
        v = v/sum(sum(D));
        if mean(v>0)>.5 || (mean(v>0)==.5 && sum(v(v>0))>-sum(v(v<0)))
            v = -v;
        end
        vectors(offset+t,:) = v';
    end
    offset = offset+lengths(s);
end
subjects = input.subjects; regions = input.regions; tr = input.tr;
input_signature = input.input_signature;
runtime_version = version;
save(output,'vectors','lengths','subjects','regions','tr','input_signature', ...
    'runtime_engine','runtime_version','-v7');
end
