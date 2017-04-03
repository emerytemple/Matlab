function [] = setAxes(S)              

% Initialization Tasks (x-direction)
range = get(S.ax,'XLim');
set(S.ed(18),'String',range(1));
set(S.ed(19),'String',range(2));

auto_ch = get(S.ax,'XLimMode');
if(strcmp(auto_ch,'auto'))
    set(S.ch(1),'Value',1);
elseif(strcmp(auto_ch,'manual'))
    set(S.ch(1),'Value',0);
end

reverse = get(S.ax,'XDir');
if(strcmp(reverse,'reverse'))
    set(S.ch(2),'Value',1);
elseif(strcmp(reverse,'normal'))
    set(S.ch(2),'Value',0);
end
    
% Initialization Tasks (y-direction)
range = get(S.ax,'YLim');
set(S.ed(20),'String',-0.5);
set(S.ed(21),'String',0.5);
set(S.ax,'YLim',[-0.5 0.5]);

auto_ch = get(S.ax,'YLimMode');
if(strcmp(auto_ch,'auto'))
    set(S.ch(3),'Value',1);
elseif(strcmp(auto_ch,'manual'))
    set(S.ch(3),'Value',0);
end

reverse = get(S.ax,'YDir');
if(strcmp(reverse,'reverse'))
    set(S.ch(4),'Value',1);
elseif(strcmp(reverse,'normal'))
    set(S.ch(4),'Value',0);
end

end