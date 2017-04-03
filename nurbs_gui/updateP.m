function [P] = updateP(S)

P.x1 = str2double(get(S.ed(1),'String'));
P.x2 = str2double(get(S.ed(3),'String'));
P.x3 = str2double(get(S.ed(5),'String'));
P.x4 = str2double(get(S.ed(7),'String'));
P.x5 = str2double(get(S.ed(9),'String'));

P.n = str2double(get(S.ed(11),'String'));
P.np = str2double(get(S.ed(12),'String'));
P.param = str2double(get(S.ed(13),'String'));
P.selectbc = str2double(get(S.ed(14),'String'));
P.clamped = str2double(get(S.ed(15),'String'));
P.tol = str2double(get(S.ed(16),'String'));
P.iter = str2double(get(S.ed(17),'String'));

end