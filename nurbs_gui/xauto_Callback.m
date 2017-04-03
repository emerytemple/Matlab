function [] = xauto_Callback(source,eventdata,S)

val = get(S.ch(1),'Value');
if(val == 0)
    set(S.ax,'XLimMode','manual');
elseif(val == 1)
    set(S.ax,'XLimMode','auto');
    
    range = get(S.ax,'XLim');
    set(S.ed(18),'String',num2str(range(1)));
    set(S.ed(19),'String',num2str(range(2)));
end

end