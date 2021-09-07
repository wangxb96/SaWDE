function configureFigure(hf,ha) 
% 作用:按论文要求配置图像各个对象的属性,如图像大小,字体,字号,颜色等

% set(hf,'Renderer','painters','PaperSize',[20.98 29.68],...
%     'InvertHardcopy','off',...
%     'Color',[1 1 1]);
% 
% set(ha,'Parent',hf,'YScale','log','YMinorTick','on','FontSize',9,...
%     'FontName','times');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 3e+005]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 6e+010]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[-1 1]);
% box(ha,'on');
% hold(ha,'all');
% grid off;

% Create figure
% set(hf,'PaperSize',[20.98/2 29.68/2]);
% oldUnit = get(hf,'unit');
% set(hf,'unit','centimeters');
% pos = get(hf,'position');
% % 首先要明白通栏和半栏的概念。一般杂志为Ａ４纸大小左右，通栏为
% % １５厘米宽，半栏为７.５厘米宽（不同杂志可能稍有不同）。因此
% % 在保存图片的时候自己就要 判断该图片是准备为通栏还是半栏发表
% % 。当然我们大部分图片是半栏（实际上杂志社也鼓励为半栏）发表，
% % 因此图片大小就在宽度设定的时候设置为7.5厘米。如 果你保存的分
% % 辨率越高，图片大小还可以更小。例如图片为半栏，分辨率为600dpi
% % ，那么宽度就可以设为4厘米。因为某图片分辩率为600dpi，它现在 
% % 的尺寸就可放大至一倍以上使用也没有问题。
% % 此处设置为宽4cm,高3cm
% pos(3:4) = [4,3];
% set(hf,'position',pos);
% set(hf,'unit',oldUnit);

% ******************************************************************
% from http://www.ymlib.net/article/sort010/info-852.html
%这4句是将字体大小改为8号字，在小图里很清晰
figure_FontSize=8;
% set(get(ha,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(ha,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',12),'FontSize',figure_FontSize);
%这句是设置绘图的大小，不需要到word里再调整大小。我给的参数，图的大小是7cm
% set(hf,'paperposition',get(hf,'position')*0.25);
% oldUnit = get(hf,'unit');
% set(hf,'unit','centimeters');
% pos = get(hf,'position');
% % 首先要明白通栏和半栏的概念。一般杂志为Ａ４纸大小左右，通栏为
% % １５厘米宽，半栏为７.５厘米宽（不同杂志可能稍有不同）。因此
% % 在保存图片的时候自己就要 判断该图片是准备为通栏还是半栏发表
% % 。当然我们大部分图片是半栏（实际上杂志社也鼓励为半栏）发表，
% % 因此图片大小就在宽度设定的时候设置为7.5厘米。如 果你保存的分
% % 辨率越高，图片大小还可以更小。例如图片为半栏，分辨率为600dpi
% % ，那么宽度就可以设为4厘米。因为某图片分辩率为600dpi，它现在 
% % 的尺寸就可放大至一倍以上使用也没有问题。
% % 此处设置为宽4cm,高3cm
% pos(3:4) = [7,5];
% set(hf,'position',pos);
% set(hf,'unit',oldUnit);
oldPos = get(hf,'Position');
oldPos(3:4) = [260,220];
set(hf,'Position',oldPos);
%这句是设置xy轴在图片中占的比例，可能需要自己微调。
% set(ha,'Position',[.13 .17 .80 .74]);
%这句是将线宽改为2
lw = 0.25;
set(findobj(get(ha,'Children'),'LineWidth',0.5),'LineWidth',lw);
%设置图片的字体类型和字号大小的。
set(ha, 'Fontname', 'Times', 'Fontsize', figure_FontSize);
% ******************************************************************

