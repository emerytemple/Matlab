function [xc, yc, xr, yr] = axisymmetric_internal(tht, xpe, rrre, gamma, nei, ne)
    me = Supersonic('Mar','x',xpe,'g',gamma);
    ve = Supersonic('nu','M',me,'g',gamma);

    vei = (ve-tht)/2.0;
    
    mei = Supersonic('Mnu','v',vei,'g',gamma);
    uei = Supersonic('mu','M',mei)*(pi/180.0);
    xpei = Supersonic('A/A*','M',mei,'g',gamma);

    ve = ve*(pi/180.0);
    vei = vei*(pi/180);
    aei = (uei+vei-ve);

    % find point P
    rpre = sqrt(1.0-mei*xpei*sin(aei)/xpe);
    xpre = (rpre-1.0)*cos(aei)/sin(aei);

    tht = pi/2.0-(tht*(pi/180.0));
    
    % calculate internal contours
    dm = (mei-1.0)/(nei-1.0);
    for i = 1:nei
        mx = 1.0+((i-1)*dm);

        ux = Supersonic('mu','M',mx)*(pi/180.0);
        vx = Supersonic('nu','M',mx,'g',gamma)*(pi/180.0);
        xpx = Supersonic('A/A*','M',mx,'g',gamma);

        bx = tht-(pi/2)-vx+abs(vei-ve);
        c = 2.0*rrre*sin(0.5*bx);
        psi = pi-tht+vx-0.5*(pi-bx);

        ax = 2.0*vei-ve-vx+ux;

        yr(i) = rpre+c*sin(psi);
        xr(i) = xpre-c*cos(psi);
        
        yc(i) = sqrt((yr(i)^2)+(mx*xpx*sin(ax)/xpe));
        xc(i) = xr(i)+(yc(i)-yr(i))*cos(ax)/sin(ax);
    end

    % calculate external contour
    dm = (me-mei)/(nei-1.0);
    for i = 1:ne
        mx = mei+((i-1)*dm);

        ux = Supersonic('mu','M',mx)*(pi/180.0);
        vx = Supersonic('nu','M',mx,'g',gamma)*(pi/180.0);
        xpx = Supersonic('A/A*','M',mx,'g',gamma);
        
        yr(nei+i) = sqrt(1-mx*xpx*sin(ve-vx+ux)/xpe);
        xr(nei+i) = (1-yr(nei+i))/tan(ve-vx+ux);
    end

    % redimensionalize
%     re = 11.864303702576;
%     for i = 1:nei
%         xc(i) = xc(i)*re;
%         yc(i) = yc(i)*re;
% 
%         xr(i) = xr(i)*re;
%         yr(i) = yr(i)*re;
%     end
%     for i = 1:ne
%         xr(nei+i) = xr(nei+i)*re;
%         yr(nei+i) = yr(nei+i)*re;
%     end
    xc'
    yc'
    xr'
    yr'
    
    plot(xr,yr,xc,yc)
end