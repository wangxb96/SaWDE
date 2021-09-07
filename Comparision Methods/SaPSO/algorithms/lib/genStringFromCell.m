function aStr = genStringFromCell(c)
% 作用:连接cell内所有字符生成字符串
% 用法: 
% a = {'a','b','c'};
% genStringFromCell(a)
% 
% ans =
% 
% a,b,c,

aStr = '';
len = length(c);
for i=1:len    
    aStr = strcat(aStr,c{i});
    aStr = strcat(aStr,', ');
end