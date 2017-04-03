function [] = xlim_Callback(source,eventdata,S)

xmin = str2double(get(S.ed(18),'String'));
xmax = str2double(get(S.ed(19),'String'));
set(S.ax,'XLim',[xmin xmax]);
set(S.ax,'XLimMode','manual');
set(S.ch(1),'Value',0);

end