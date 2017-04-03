function writetoNCfile(fid,cid,k,v,d)
    fprintf(fid,'%d\n',cid);
    fprintf(fid,'%d\n',k);
    fprintf(fid,' %d\n',length(d)-1);
    for i = 1:length(v)
        fprintf(fid,'%f\n',v(i));
    end
    for i = 1:length(d)
        fprintf(fid,'%f %f %f %f\n',d(i,1),d(i,2),0,1);
    end
end