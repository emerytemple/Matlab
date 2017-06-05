function [xc, yc, xr, yr] = first_pts(tht, xpe, rrre, gamma)
    format long

    me = Supersonic('Mar','x',xpe,'g',gamma);
    ve = Supersonic('nu','M',me,'g',gamma);

    vei = (ve-tht)/2.0;
    
    mei = Supersonic('Mnu','v',vei,'g',gamma);
    uei = Supersonic('mu','M',mei)*(pi/180.0);
    xpei = Supersonic('A/A*','M',mei,'g',gamma);

    ve = ve*(pi/180.0);
    vei = vei*(pi/180);
    aei = (uei+vei-ve);

    rpre = sqrt(1.0-mei*xpei*sin(aei)/xpe);
    xpre = (rpre-1.0)*cos(aei)/sin(aei);

    tht = pi/2.0-(tht*(pi/180.0));

    bx = tht-(pi/2)+abs(vei-ve);
    c = 2.0*rrre*sin(0.5*bx);
    psi = pi-tht-0.5*(pi-bx);

    ax = 2.0*vei-ve+(pi/2.0);

    yr = rpre+c*sin(psi);
    xr = xpre-c*cos(psi);

    yc = sqrt((yr^2)+(sin(ax)/xpe));
    xc = xr+(yc-yr)*cos(ax)/sin(ax);
end