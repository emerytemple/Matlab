function [] = keepLT(S)

% keep x1 < x5
x1 = str2double(get(S.ed(1),'String'));
x5 = str2double(get(S.ed(9),'String'));

if(x1 > x5)
    tmp = x1;
    x1 = x5;
    x5 = tmp;
end

% update edit boxes
set(S.ed(1),'String',num2str(x1));
set(S.ed(9),'String',num2str(x5));

% update sliders
set(S.sp(1),'Value',x1);
set(S.sp(5),'Value',x5);

end