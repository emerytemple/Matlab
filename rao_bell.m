function [x,y] = rao_bell(me, gamma, np)
    ae = Supersonic('mu','M',me)*(pi/180.0)
    aed = ae*(180.0/pi)
    the = asin((2.0/(gamma*(me^2)))*cot(ae))/2.0
    thed = the*(180.0/pi)
    mes = sqrt(1.0/(gamma-1.0+(2.0/(me^2))))
    
    ppoe = 1/Supersonic('Po/P','M',me,'g',gamma);
    bot = (me^2)*ppoe*(sin(the)^2)*tan(ae);
    
    dm = (me-1.103)/(np-1.0);
    for i = 1:np
        mx = 1.103+((i-1)*dm);
        
        ms = sqrt(1.0/(gamma-1.0+(2.0/(mx^2))));
        ax = Supersonic('mu','M',mx)*(pi/180.0);
        ad = ax*(180.0/pi);
        
%         fth = @(th) ((ms*(cos(th-a)/cos(a)))-(mes*(cos(the-ae)/cos(ae))))
%         th = fzero(fth, 1)
%         th = acos((mes/ms)*(cos(a)/cos(ae))*cos(the-ae))+ a
%         thd = th*(180.0/pi)
        
        thx = (-(6.966*(mx^2))+(23.687*mx)+15.203)*(pi/180.0);
        
        ppox = 1/Supersonic('Po/P','M',mx,'g',gamma);
        top = (mx^2)*ppox*(sin(thx)^2)*tan(ax);

        yye = top/bot;
        
        ansss(i,1) = mx;
        ansss(i,2) = yye;
        
    end
    ansss
end

