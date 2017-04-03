function [] = ylim_Callback(source,eventdata,S)

ymin = str2double(get(S.ed(20),'String'));
ymax = str2double(get(S.ed(21),'String'));
set(S.ax,'YLim',[ymin ymax]);
set(S.ax,'YLimMode','manual');
set(S.ch(3),'Value',0);

end