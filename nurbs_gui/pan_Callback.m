function [] = pan_Callback(source,eventdata,S)

if(get(S.tg(2),'Value') == 1)
    set(S.tg(2),'Value',0);    
end

pan;

setAxes(S);

end