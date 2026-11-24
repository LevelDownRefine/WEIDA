function filtered = original_preprocess(BOLD, TR)
% Original WEIDA-DBS amplitude preprocessing, wrapped for numerical testing.
% Source: https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L88
% Original statements: lines 70-76 and 87-88; algorithm statements unchanged.
% Adaptation: function arguments, region loop and amplitude output only.
% Hilbert phase extraction and all subsequent WEIDA steps are excluded.
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k
filtered = zeros(size(BOLD));
for seed=1:size(BOLD,1)
    ts=demean(detrend(BOLD(seed,:)));
        signal_filt =filtfilt(bfilt,afilt,ts);
    filtered(seed,:) = signal_filt;
end
end
