function aStr = genStringFromVec(v)
% 作用:连接cell内所有字符生成字符串
% 用法: 
% a = [1,2,3];
% genStringFromVec(a)
% 
% ans =
% 
% 1,2,3,

aStr = '';
len = length(v);
for i=1:len    
    aStr = strcat(aStr,num2str(v(i)));
    aStr = strcat(aStr,', ');
end