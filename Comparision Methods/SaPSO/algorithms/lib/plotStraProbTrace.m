function plotStraProbTrace(fh,xlsFileName) 
% 作用:画收敛图,输入为xls文件,其数据组织格式为:第一列为FEs,第二列为对应的erro value或者最优值,第一第二列为对应算法1,以此类推

x = xlsread(xlsFileName);

[r,l] = size(x);
totalAlgorithm = l/2;
totalStrategy = 4;

ph = {};

styleLine = {'-','--',':','-.'};
styleMaker = {'x','o','*','d','+','.','none','p','s','<','>'}; % 可选的点的样式
styleColor = 'k'; % balck:k {'k','r','g'} 曲线颜色
algorithm = {'CLONALG','opt-IA','BCSA','CLPSO','saDE','CMA-ES','GL-25','SALIA'};  % 算法名称,按序对应于每条曲线的legend说明
strategy = {'Mr1','Mcr1','ISR','Mr2'}


oldYmin = [];
oldYmax = [];
for i=1:totalStrategy
    ph{i} = plot(x(:,1),x(:,i+1));hold on; 
    set(ph{i},'LineStyle', '-','Color',styleColor,'Marker',styleMaker{i},'DisplayName',strategy{i});
    v = axis;
    oldYmin = [oldYmin v(3)];
    oldYmax = [oldYmax v(4)];
end

% 设定X轴,Y轴的标题
ylabel('Probability');
xlabel('FEs');
% 设定X轴,Y轴的刻度范围

axis([0,300000,min(oldYmin),max(oldYmax)]);
% set(gca,'YScale','log');
%set(gca,'Xtick',[0:50000:300000]);

% 显示各图线的legend
legend('show')

fileName = strcat('a_pic',xlsFileName(1:end-4));
saveas(gcf,fileName,'eps')
saveas(gcf,fileName,'tif')
saveas(gcf,fileName,'fig')
