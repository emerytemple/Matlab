function [] = read_Callback(source,eventdata,S,P)

[filename,pathname] = uigetfile('*.dat','File Selector');

% read in airfoil data
fulll = fullfile(pathname,filename);
fid1 = fopen(fulll,'r');
points = fscanf(fid1, '%g %g', [2 inf]);
status = fclose(fid1);
myData.pts = points';

myData.k = 4;
myData.clamped = 1;

guidata(S.fh,myData)

plotS809(S);

end