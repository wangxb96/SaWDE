function strategySeq = getStrategyOnAParticle(P,num)
% 作用:轮盘赌方式选择进化策略,num为抛筛子次数 
% 参考: http://blog.163.com/coldstar@126/blog/static/2317011200681792853432/
% 修改时间: 2010-11-16 
% 改动:
% (1)增加确定选择步骤
% -----------------------------------------------
% 使用示例
% strategySeq = getStrategyOnAParticle([0.25,0.25,0.25,0.25],4*10)
% strategySeq =
%     2
% **********************************************************

% 精度太低
% m = length(P); % 若P为4维向量,分量和为1,m=4
% Select = zeros(1,num);
% for i=1:num
%     Select(i) = m;% 初始化为最后一个
%     for j=1:m %:按概率选择
%         if P(j)>rand()
%             Select(i)=j;
%             break;
%         end
%     end
% end
% 精度太低


m = length(P); % 若P为4维向量,分量和为1,m=4
Select = zeros(1,num); % Select = [0     0     0     0     0]
r = rand(1,num); % r = [0.6443    0.3786    0.8116    0.5328    0.5328]
for i=1:num
    sumP = 0;
    j = ceil(m*rand); %产生1~m之间的随机整数
    while sumP < r(i)
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end
    %Select(i) = mod(j-1,m)+1-1;
    Select(i) = mod(j-2,m)+1;
end

% 确定最终选择
tempSelect = zeros(1,m);
for i=1:m
    tempSelect(i) = sum(Select==i);
end
[val,ind] = sort(tempSelect,'descend');
strategySeq = ind(1);