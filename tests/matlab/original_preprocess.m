function filtered = original_preprocess(BOLD, TR)
% 将原版 WEIDA-DBS 幅值预处理封装为函数，用于数值对照测试。
% 出处： https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L88
% 原语句位于第 70–76、87–88 行；算法语句保持不变。
% 仅添加函数参数、脑区循环及幅值输出；注释翻译为中文。
% 不包含 Hilbert 相位提取及后续 WEiDA 步骤。
fnq=1/(2*TR);                 % 奈奎斯特频率
flp = 0.04;                    % 带通下限频率（Hz）
fhi = 0.07;                    % 带通上限频率（Hz）
Wn=[flp/fnq fhi/fnq];         % Butterworth 带通的归一化频率
k=2;                          % Butterworth 滤波器的阶数参数为 2
[bfilt,afilt]=butter(k,Wn);   % 构建滤波器
clear fnq flp fhi Wn k
filtered = zeros(size(BOLD));
for seed=1:size(BOLD,1)
    ts=demean(detrend(BOLD(seed,:)));
        signal_filt =filtfilt(bfilt,afilt,ts);
    filtered(seed,:) = signal_filt;
end
end
