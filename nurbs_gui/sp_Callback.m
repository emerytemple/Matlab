function [] = sp_Callback(source,eventdata,S)

val = get(source,'Value');

for i = 1:length(S.sp)
    if(source == S.sp(i))
        ind = (2*i)-1;
        break;
    end
end

set(S.ed(ind),'String',num2str(val));

keepLT(S);

plotS809(S);

end