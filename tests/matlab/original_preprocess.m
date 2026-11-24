function filtered = original_preprocess(BOLD, TR)
% Test-only wrapper; execute unchanged preprocessing from the upstream script.
% Source and Git blob hashes: upstream/WEIDA-DBS/SOURCE.md.
upstream = fullfile(fileparts(mfilename('fullpath')), 'upstream', 'WEIDA-DBS');
source = fileread(fullfile(upstream, 'StateSpace_modelling.m'));
settings = regexp(source, '(?m)^fnq=1/\(2\*TR\);[\s\S]*?^clear fnq flp fhi Wn k', 'match');
amplitude = regexp(source, '(?m)^\s*ts=demean\(detrend\(BOLD\(seed,:\)\)\);\s*signal_filt =filtfilt\(bfilt,afilt,ts\);', 'match');
assert(numel(settings)==1 && numel(amplitude)==1, ...
    'Expected exactly one original filter setup and amplitude block');
oldPath = path;
cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(upstream, 'functions'), '-begin');
assert(strcmp(which('demean'), fullfile(upstream, 'functions', 'demean.m')));
eval(settings{1});
filtered = zeros(size(BOLD));
for seed=1:size(BOLD,1)
    eval(amplitude{1});
    filtered(seed,:) = signal_filt;
end
end
