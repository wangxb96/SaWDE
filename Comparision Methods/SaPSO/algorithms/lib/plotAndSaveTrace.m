% plotAndSaveTrace.m
% 
% xlsFileName = {'A_strace_WangBenchmarkFun_1.xls'
% 'A_strace_WangBenchmarkFun_2.xls'
% 'A_strace_WangBenchmarkFun_3.xls'
% 'A_strace_WangBenchmarkFun_4.xls'
% 'A_strace_WangBenchmarkFun_5.xls'
% 'A_strace_WangBenchmarkFun_6.xls'
% 'A_strace_WangBenchmarkFun_7.xls'
% 'A_strace_WangBenchmarkFun_8.xls'
% 'A_strace_WangBenchmarkFun_9.xls'
% 'A_strace_WangBenchmarkFun_10.xls'
% 'A_strace_WangBenchmarkFun_11.xls'
% 'A_strace_WangBenchmarkFun_12.xls'
% 'A_strace_WangBenchmarkFun_13.xls'
% 'A_strace_WangBenchmarkFun_14.xls'
% 'A_strace_WangBenchmarkFun_15.xls'
% 'A_strace_WangBenchmarkFun_16.xls'
% 'A_strace_WangBenchmarkFun_17.xls'
% 'A_strace_WangBenchmarkFun_18.xls'
% 'A_strace_WangBenchmarkFun_19.xls'
% 'A_strace_WangBenchmarkFun_20.xls'
% 'A_strace_WangBenchmarkFun_21.xls'
% 'A_strace_WangBenchmarkFun_22.xls'
% 'A_strace_WangBenchmarkFun_23.xls'
% 'A_strace_WangBenchmarkFun_24.xls'
% 'A_strace_WangBenchmarkFun_25.xls'
% 'A_strace_WangBenchmarkFun_26.xls'};

% for i=[1:26]
%     fh{i} = figure(i);
%     plotTrace(fh{i},xlsFileName{i}); 
% end

% xlsFileName2 = {'A_straSelProb_WangBenchmarkFun_01.xls'
% 'A_straSelProb_WangBenchmarkFun_05.xls'
% 'A_straSelProb_WangBenchmarkFun_08.xls'
% 'A_straSelProb_WangBenchmarkFun_17.xls'}
% 
xlsFileName3 = {'A_strace_WangBenchmarkFun_1.xls'
'A_strace_WangBenchmarkFun_2.xls'
'A_strace_WangBenchmarkFun_3.xls'
'A_strace_WangBenchmarkFun_4.xls'
'A_strace_WangBenchmarkFun_5.xls'
'A_strace_WangBenchmarkFun_6.xls'
'A_strace_WangBenchmarkFun_7.xls'
'A_strace_WangBenchmarkFun_8.xls'
'A_strace_WangBenchmarkFun_9.xls'
'A_strace_WangBenchmarkFun_10.xls'
'A_strace_WangBenchmarkFun_11.xls'
'A_strace_WangBenchmarkFun_12.xls'
'A_strace_WangBenchmarkFun_13.xls'
'A_strace_WangBenchmarkFun_14.xls'
'A_strace_WangBenchmarkFun_15.xls'
'A_strace_WangBenchmarkFun_16.xls'}

% for i=1:length(xlsFileName2)
%     fh{i} = figure(i);
%     plotStraProbTrace(fh{i},xlsFileName2{i}); 
% end
for i=[2,4,6,8,9,15]
    fh{i} = figure(i);
    configureFigure(fh{i},gca)
    plotTrace(fh{i},xlsFileName3{i}); 
end

% xlsFileName2 = {'A_straSelProb_WangBenchmarkFun_01.xls'
% 'A_straSelProb_WangBenchmarkFun_05.xls'
% 'A_straSelProb_WangBenchmarkFun_08.xls'
% 'A_straSelProb_WangBenchmarkFun_17.xls'}

% xlsFileName3 = {'vecNum_1.xls'
%     'vecNum_7.xls'}

% for i=[1:2]
%     fh{i} = figure(i);
%     configureFigure(fh{i},gca)
%     plotVecNum(fh{i},xlsFileName3{i}); 
% end


