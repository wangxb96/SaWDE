function plotTrace(fh,xlsFileName) 
% 作用:画收敛图,输入为xls文件,其数据组织格式为:第一列为FEs,第二列为对应的erro value或者最优值,第一第二列为对应算法1,以此类推

x = xlsread(xlsFileName);

[r,l] = size(x);
totalAlgorithm = l/2;

ph = {};

styleLine = {'-','--',':','-.'};
styleMaker = {'.','o','*','d','+','x','none','p','s','<','>'}; % 可选的点的样式
styleColor = 'k'; % balck:k {'k','r','g'} 曲线颜色
% algorithm = {'CLONALG','opt-IA','BCSA','CLPSO','saDE','CMA-ES','GL-25','SALIA'};  % 算法名称,按序对应于每条曲线的legend说明
algorithm = {'CLONALG','BCSA','MBCSA'};


oldYmin = [];
oldYmax = [];
for i=1:totalAlgorithm
    xx=x(:,(i-1)*2+1)/(10.^5);
    ph{i} = plot(xx,x(:,i*2));hold on; 
    set(ph{i},'Color',styleColor,...
        'Marker',styleMaker{i},...
        'DisplayName',algorithm{i});
    v = axis;
    oldYmin = [oldYmin v(3)];
    oldYmax = [oldYmax v(4)];
end

% 设定X轴,Y轴的标题
ylabel('log(\it{error})');
xlabel('Funtion Evaluations /10^5');
% 设定X轴,Y轴的刻度范围

axis([0,3,min(oldYmin),max(oldYmax)]);
set(gca,'YScale','log');
%set(gca,'Xtick',[0:50000:300000]);

% 显示各图线的legend
legend('show','Location','SouthWest')
legend('boxoff')

fileName = strcat('a_pic',xlsFileName(1:end-4));
% saveas(gcf,fileName,'eps')
% saveas(gcf,fileName,'tif')
saveas(gcf,fileName,'fig')
