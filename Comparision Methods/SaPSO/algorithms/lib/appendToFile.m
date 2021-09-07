function appendToFile(fileName,curIter,fitcount,gbestVal)
% 作用: 将相关结果添加到文件fileName后面

%计算平均值
fid = fopen(fileName,'a+');
fprintf(fid,'%d\t%d\t%1.20E\t\r\n',curIter,fitcount,gbestVal);
fclose(fid);