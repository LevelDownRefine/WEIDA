function run_reference()
% Run this function in MATLAB or Octave. No Python or real data required.
root=fileparts(mfilename('fullpath'));
previous=pwd; cd(root);
cleanup=onCleanup(@() cd(previous));
diary('runtime.log');
fprintf('Runtime: %s\n',version);
ver;
try
    reference;
    fprintf('SUCCESS: return reference.mat and runtime.log with input and manifest.\n');
catch err
    fprintf('FAILED: %s\n',err.message);
    diary off;
    rethrow(err);
end
diary off;
end
