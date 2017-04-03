function [] = write_Callback(source,eventdata,S)

[filename,pathname] = uiputfile('*.NC','Save Airfoil As');
fulll = fullfile(pathname,filename);

k = str2double(get(S.ed(11),'String'));

myData = guidata(S.fh);
u1 = myData.u1;
ucc = myData.ucc;
u2 = myData.u2;

d1 = myData.d1;
cc = myData.cc;
d2 = myData.d2;

fid2 = fopen(fulll,'w'); % open file for writing
writetoNCfile(fid2,1,k,u1,d1); % segment 1
writetoNCfile(fid2,2,k,ucc,cc); % cubic curve
writetoNCfile(fid2,3,k,u2,d2); % segment 2
status = fclose('all');

end