function [] = sto_Callback(source,eventdata,S)

% find which x variable is being changed
for i = 1:length(S.sp)
    if(source == S.pb(i+2))
        ind = (2*i)-1;
        ind2 = 3*i;
        break;
    end
end

xval = get(S.ed(ind),'String');
stoval = get(S.ed(ind+1),'String');
stonum = str2double(stoval);

% get bounds
minstr = get(S.tx(ind2-1),'String');
minnum = str2double(minstr);
maxstr = get(S.tx(ind2),'String');
maxnum = str2double(maxstr);

% check bounds & validation check (make a number)
if((stonum < minnum) || isnan(stonum))
    stoval = minstr;
    stonum = minnum;
elseif(stonum > maxnum)
    stoval = maxstr;
    stonum = maxnum;
end 

% set values
set(S.ed(ind),'String',stoval);
set(S.ed(ind+1),'String',xval);
set(S.sp(i),'Value',stonum);

keepLT(S);

plotS809(S);

end