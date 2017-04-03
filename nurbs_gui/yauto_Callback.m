function [] = yauto_Callback(source,eventdata,S)

val = get(S.ch(3),'Value');
if(val == 0)
    set(S.ax,'YLimMode','manual');
elseif(val == 1)
    set(S.ax,'YLimMode','auto');
    
    range = get(S.ax,'YLim');
    set(S.ed(20),'String',num2str(range(1)));
    set(S.ed(21),'String',num2str(range(2)));
end

end