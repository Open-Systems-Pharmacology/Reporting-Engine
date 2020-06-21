function printFigureHighResPng(figH,savePath)


rez=300; %resolution (dpi) of final graphic
f=figH; %f is the handle of the figure you want to export
figpos=getpixelposition(f); 
resolution=get(0,'ScreenPixelsPerInch'); 
set(f,'paperunits','inches','papersize',...
    figpos(3:4)/resolution,'paperposition',...
    [0 0 figpos(3:4)/resolution]); 
set(f,'PaperOrientation','portrait');
print(f,fullfile(savePath),'-dpng',['-r',num2str(rez)],'-opengl') %save file