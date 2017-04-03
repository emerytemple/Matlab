function [] = zoom_Callback(source,eventdata,S)

if(get(S.tg(1),'Value') == 1)
    set(S.tg(1),'Value',0);
end

zoom;

setAxes(S);

end