function BOLD = original_read_bold(subpth)
% 出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L25-L28
% 保留作者读取语句，仅增加文件路径参数并将 ts_cell 输出命名为 BOLD。
a = importdata(subpth);
a = a.data;
a(:,1)=[];
BOLD = a';
end
