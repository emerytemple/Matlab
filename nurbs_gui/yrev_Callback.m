function [] = yrev_Callback(source,eventdata,S)

val = get(S.ch(4),'Value');
if(val == 0)
    set(S.ax,'YDir','normal');
elseif(val == 1)
    set(S.ax,'YDir','reverse');
end

end