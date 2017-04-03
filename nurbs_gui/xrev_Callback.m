function [] = xrev_Callback(source,eventdata,S)

val = get(S.ch(2),'Value');
if(val == 0)
    set(S.ax,'XDir','normal');
elseif(val == 1)
    set(S.ax,'XDir','reverse');
end

end