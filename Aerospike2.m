function Aerospike2()
wtf1 = fzero(@(x) root2(x), 15);
[xc1, yc1, xr1, yr1] = axisymmetric_internal(2,15,90,wtf1,1.4,100,100);
check_linear_internal(xc1, yc1, xr1, yr1)
end

function [ar] = root2(crapval)
    [xc,yc,xr,yr] = axisymmetric_internal(2,15,90,crapval,1.4,100,100);

    re = abs(yr(length(yr))) - abs(yc(length(yc)));
    Ae = pi*(re^2);
    
    re = abs(yr(length(yr)));
    r = re - abs(yr(1));
    alpha = (pi/2)-atan((abs(xr(1)))/(abs(yr(1))));
    At = abs((pi/sin(alpha))*((re^2)-(r^2)));

    throat_len = abs(re-r)/sin(alpha);
    throat_len2 = sqrt(((yc(1)-yr(1))^2)+((xc(1)-xr(1))^2));
    ar = Ae/At-175
end

function [xc, yc, xr, yr, xp, yp] = axisymmetric_internal(xpei, xpe, tht, rrre, gamma, nei, ne)
    tht = tht*(pi/180.0);

    me = Supersonic('Mar','x',xpe,'g',gamma);
    ve = Supersonic('nu','M',me,'g',gamma)*(pi/180.0);
    
    mei = Supersonic('Mar','x',xpei,'g',gamma);
    vei = Supersonic('nu','M',mei,'g',gamma)*(pi/180.0);
    uei = Supersonic('mu','M',mei)*(pi/180.0);
    xpei = Supersonic('A/A*','M',mei,'g',gamma);

    aei = (uei+vei-ve);

    % find point P
    rpre = sqrt(1.0-mei*xpei*sin(aei)/xpe);
    xpre = (rpre-1.0)*cos(aei)/sin(aei);
    
    xp = xpre;
    yp = rpre;
    
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
    rxre(1) = sqrt(1-(mei*xpei)*sin(ve-vei+uei)/xpe); % fix this later
    xxre(1) = (1-rxre(1))/tan(ve-vei+uei);

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
    re = 11.864303702576;
    for i = 1:nei
        xc(i) = xc(i)*re;
        yc(i) = yc(i)*re;

        xr(i) = xr(i)*re;
        yr(i) = yr(i)*re;
    end
    for i = 1:ne
        xr(nei+i) = xr(nei+i)*re;
        yr(nei+i) = yr(nei+i)*re;
    end
end




























% variables:
% m - mach number
% v - prandtl-meyer function
% u - mach angle
% a - angle flow makes with horizontal
% xp - area ratio
% n - number of points
% tht - throat angle
%
% subscripts:
% e - end of external
% ei - end of internal
% x - current