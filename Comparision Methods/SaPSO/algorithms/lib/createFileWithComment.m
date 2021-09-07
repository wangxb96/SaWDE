function fid = createFileWithComment(fileName,comment)
% 作用: 将相关结果添加到文件fileName后面

%计算平均值
fid = fopen(fileName,'a+');
fprintf(fid,'%s',comment);
