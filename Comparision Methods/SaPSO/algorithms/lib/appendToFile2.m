function appendToFile2(fileName,gbestMean,gbestStd,fitcountMean,fitcountStd)
% 作用: 将相关结果添加到文件fileName后面

%计算平均值
fid = fopen(fileName,'a+');
fprintf(fid,'%E\t%E\t%E\t%E\t\r\n',gbestMean,gbestStd,fitcountMean,fitcountStd);
fclose(fid);