function[] = Aerospike()
clear all
clc
format long

% [x1,y1] = linear_external(20,1,1.4,20);
% check_linear_external(x1,y1)

% [x2,y2] = axisymmetric_external(20,1,1.4,20);
% check_axisymmetric_external(x2,y2)

% [x3,y3] = axisymmetric_external2(20,1,1.4,20);
% check_axisymmetric_external(x3,y3)

% plot(x1,y1,x2,y2,x3,y3); %,x4,y4);

% wtf1 = fzero(@(x) root1(x), 1);
% [xc1, yc1, xr1, yr1] = linear_internal(200,1,90,1.2575,100,100,wtf1) % here
% check_linear_internal(xc1, yc1, xr1, yr1)

wtf2 = fzero(@(x) root2(x), 175);
[xc2, yc2, xr2, yr2] = axisymmetric_internal(wtf2,1,90,1.2575,50,75,0.17);
check_axisymmetric_internal(xc2, yc2, xr2, yr2)

% hold on
% plot(xc1,yc1,xr1,yr1)
plot(xc2,yc2,xr2,yr2)
% hold off


end

function [ar] = root1(crapval)
    [xc, yc, xr, yr] = linear_internal(200,1,90,1.2575,100,100,crapval); % here

    throat_len = sqrt(((yc(1)-yr(1))^2)+((xc(1)-xr(1))^2));
    re = abs(yr(length(yr))) - abs(yc(length(yc)));
    ar = re/throat_len-200 % and finally here
end

function [ar] = root2(crapval)
    [xc,yc,xr,yr] = axisymmetric_internal(crapval,1,90,1.2575,100,100,0.17);

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

% linear contours
function [x,y] = linear_external(xp, tlen, gamma, nr) % linear external (approximate) contour
% Angelino, G. - Approximate Method for Plug Nozzle Design (AIAA)

    x = zeros(nr,1);
    y = zeros(nr,1);
    
    mach_exit = Supersonic('Mar','x',xp,'g',gamma);
    theta = Supersonic('nu','M',mach_exit,'g',gamma); % ve
    dm = (mach_exit-1.0)/(nr-1.0);
    
    for i = 1:nr
        M = 1.0+((i-1)*dm);
        mu = Supersonic('mu','M',M);
        nu = Supersonic('nu','M',M,'g',gamma);
        alpha = (mu+theta-nu)*(pi/180.0);
        xsi = M*Supersonic('A/A*','M',M,'g',gamma);

        % redimensionalize
        L = xsi*tlen;
        x(i) = L*cos(alpha);
        y(i) = -L*sin(alpha);
    end
end

function [xc, yc, xr, yr] = linear_internal(xp, tlen, tang, gamma, n1, n2, beta)

Me = Supersonic('Mar','x',xp,'g',gamma);
ve = Supersonic('nu','M',Me,'g',gamma);

tht = 90 - tang;
vei = 0.5*(tht+ve);
% Mei = Supersonic('Mnu','v',vei,'g',gamma);

% Mei = Supersonic('Mar','x',xpei,'g',gamma);
% vei = Supersonic('nu','M',Mei,'g',gamma);
% tht = ve - (2*vei);

cowl(1,1) = 0.0;
cowl(1,2) = 0.0;

ramp(1,1) = -tlen*sind(abs(tht));
ramp(1,2) = tlen*cosd(tht);

P(1) = Supersonic('Po/P','M',1.0,'g',gamma);

dv = vei/(n1-1);
for i = 2:n1
    v = (i-1)*dv;
    th = tht + ((i-1)*dv); % dv = dth
    
    M = Supersonic('Mnu','v',v,'g',gamma);
    P(i) = Supersonic('Po/P','M',M,'g',gamma);
    alpha = Supersonic('mu','M',M);

    angle = abs(tht) + (v/2);
    len = sqrt(2 * ((beta*tlen)^2) * (1-cosd(v)));
    ramp(i,1) = ramp(1,1) + (len * cosd(angle));
    ramp(i,2) = ramp(1,2) + (len * sind(angle));
    
    A = [1, -tand(th-dv); 1, -tand(th-alpha)];
    B = [cowl(i-1,2) - tand(th-dv)*cowl(i-1,1);
        ramp(i,2) - tand(th-alpha)*ramp(i,1)];
    solution = A\B;
    cowl(i,2) = solution(1,1);
    cowl(i,1) = solution(2,1);
end

dv = (ve-vei)/(n2-1);
th = th + dv;
ind = length(cowl);

% ls = [ramp(ind,1), ramp(ind,2); cowl(ind,1), cowl(ind,2)];
% plot(ls(:,1),ls(:,2))

for i = 1:n2
    v = vei + ((i-1)*dv);
    th = th - dv;

    M = Supersonic('Mnu','v',v,'g',gamma);
    P(ind+i) = Supersonic('Po/P','M',M,'g',gamma);
    alpha = Supersonic('mu','M',M);

    A = [1 -tand(th); 1 -tand(th+alpha)];
    B = [ramp((ind+i-1),2) - tand(th) * ramp((ind+i-1),1);
        cowl(ind,2) - (tand(th+alpha) * cowl(ind,1))];

    solution = A\B;
    ramp(ind+i,2) = solution(1,1);
    ramp(ind+i,1) = solution(2,1);
end

len = ramp(length(ramp),2);
ramp(:,2) = len - ramp(:,2);
cowl(:,2) = len - cowl(:,2);

xc = cowl(:,1);
yc = cowl(:,2);
xr = ramp(:,1);
yr = ramp(:,2);

    dist = yc(1);
    for i = 1:length(xc)
        yc(i) = yc(i) - dist;
    end
    for i = 1:length(xr)
        yr(i) = yr(i) - dist;
    end
end

% axi-symmetric contours
function [x,y] = axisymmetric_external(xp, tlen, gamma, nr) % axi-symmetric external contour
% angelino

    x = zeros(nr,1);
    y = zeros(nr,1);

    mach_exit = Supersonic('Mar','x',xp,'g',gamma);
    dm = (mach_exit-1.0)/(nr-1.0);

    theta = Supersonic('nu','M',mach_exit,'g',gamma); % ve
    theta_rad = theta*(pi/180.0);

    re = xp*tlen*(1+sqrt(1-(cos(theta_rad)/xp)));
    
    for i = 1:nr
        M = 1.0+((i-1)*dm);
        mu = Supersonic('mu','M',M);
        nu = Supersonic('nu','M',M,'g',gamma);
        alpha = (mu+theta-nu)*(pi/180.0);
        eps = Supersonic('A/A*','M',M,'g',gamma);
        xsi = (1.0-sqrt(1.0-((eps/xp)*M*sin(alpha))))/sin(alpha);

        % redimensionalize
        L = xsi*re;
        x(i) = L*cos(alpha);
        y(i) = -L*sin(alpha);
    end
end

function [x, y] = axisymmetric_external2(xp, tlen, gama, nr) % axi-symmetric external contour 2
% lee & thompson

    xm = zeros(nr,1);
    rxre = zeros(nr,1);
    xxre = zeros(nr,1);

    rm = Supersonic('Mar','x',xp,'g',gama);
    ve = Supersonic('nu','M',rm,'g',gama)*(pi/180.0);

    theta = (pi/2.0)-ve;
    ht = (xp-sqrt(xp*(xp-sin(theta))))/(xp*sin(theta));

    drm = (rm-1.0)/nr;

    xm(1)=1.0;
    rxre(1) = 1.0-ht*sin(theta);
    xxre(1) = -ht*cos(theta);

    i = 2;
    while i-nr-1 <= 0
        xm(i) = xm(i-1) + drm;
        vx = Supersonic('nu','M',xm(i),'g',gama)*(pi/180.0);
        y = 1.0/xm(i);
        ux = atan(y/sqrt(1.0-y*y));

        a2 = (gama+1.0)/(2.0*(gama-1.0));
        ToT = Supersonic('To/T','M',xm(i),'g',gama);
        b2 = (2.0/(gama+1.0))*ToT;
        rxre(i) = sqrt(1.0-(b2^a2)*sin(ve-vx+ux)/xp);
        xxre(i) = (1.0-rxre(i))*cos(ve-vx+ux)/sin(ve-vx+ux);

        i = i+1;
    end

    x = real(xxre)';
    y = real(rxre)'-1.0;
    
    alen = sqrt((x(1)^2)+(y(1)^2));
    scale = tlen/alen;
    
    x = scale.*x;
    y = scale.*y;
end

function [xc, yc, xr, yr] = axisymmetric_internal(xp, tlen, tang, gamma, n1, n2, rrre) % axi-symmetric internal-external contour
% lee & thompson
% n1 - number of discretization points for internal section
% n2 - number of discretization points for external section
% xp - end of external area ratio
% xp2 - end of internal area ratio
% rrre - nondimensional radius of circle defining the throat
% pht - angle (in degrees) the throat makes with horizontal
% gamma - ratio of specific heats

    tang = tang*(pi/180.0);

    gp1 = gamma+1.0;
    gm1 = gamma-1.0;
    exp = gp1/(2.0*gm1);
    c1 = 2.0/gp1;
    c3 = gp1/gm1;
    c4 = (1.0-gamma)/gamma;
    c5 = 2.0/gm1;

    rm = Supersonic('Mar','x',xp,'g',gamma);
    ve = Supersonic('nu','M',rm,'g',gamma)*(pi/180.0);
    
    % rmei = sqrt(c5*((peipc^c4)-1.0));
    rmei = Supersonic('Mar','x',2.0,'g',gamma);
    P(1) = Supersonic('Po/P','M',rmei,'g',gamma);
    vei = Supersonic('nu','M',rmei,'g',gamma)*(pi/180.0);
%     tht = 90 - tang;
%     vei = 0.5*(tht+ve);
%     rmei = Supersonic('Mnu','v',vei,'g',gamma);

    thei = vei - ve;
    phei = thei + atan((1/rmei)/sqrt(1-((1/rmei)*(1/rmei))));

    rpre = sqrt(1-(c1*(1+(0.5*gm1*rmei*rmei)))*sin(phei)/xp);
    xpre = (rpre-1)*cos(phei)/sin(phei);
    drm = (rmei-1.0)/n1;

    % internal portion of the nozzle (1 is ramp, 2 is cowl)
    for k = 1:n1+1
        if k == 1
            xm(1) = 1;
        end

        vx = Supersonic('nu','M',xm(k),'g',gamma)*(pi/180.0);
        bx = tang-(pi/2)-vx+abs(thei);
        x1re = 2.0*rrre*sin(0.5*bx);
        psi = pi-tang+vx-0.5*(pi-bx);
        rx1re(k) = rpre+x1re*sin(psi);
        xx1re(k) = xpre-x1re*cos(psi);

        if k == 1
            rx2re(1) = sqrt(rx1re(1)*rx1re(1)+sin(tang)/xp);
            xx2re(1) = xx1re(1)+(rx2re(1)-rx1re(1))*cos(tang)/sin(tang); % use atan
        else
            phx = 2.0*vei-ve-vx+atan(1.0/(xm(k)*sqrt(1.0-(1.0/xm(k))^2.0)));
            a2 = c1*(1.0+0.5*gm1*xm(k)*xm(k));
            rx2re(k) = sqrt(rx1re(k)*rx1re(k)+(a2^(0.5*c3))*sin(phx)/xp);
            xx2re(k) = xx1re(k)+(rx2re(k)-rx1re(k))*cos(phx)/sin(phx);
        end
        if k < n1+1
            xm(k+1) = xm(k) + drm;
            P(k+1) = Supersonic('Po/P','M',xm(k+1),'g',gamma);
        end
        
    end

%     ind = length(xx1re);
%     ls = [xx1re(ind), rx1re(ind); xx2re(ind), rx2re(ind)];
%     plot(ls(:,1),ls(:,2))

    % external portion of the nozzle
    xm2(1) = xm(k);

    ux = atan(1.0/(xm2(1)*sqrt(1.0-(1.0/xm2(1))^2.0)));
    vx = Supersonic('nu','M',xm2(1),'g',gamma)*(pi/180.0);
    rxre(1) = sqrt(1-((c1*(1+0.5*gm1*xm2(1)*xm2(1)))^exp)*sin(ve-vx+ux)/xp);
    xxre(1) = (1-rxre(1))*cos(ve-vx+ux)/sin(ve-vx+ux);

    drm = (rm-xm2(1))/n2;
    xm2(2) = xm2(1)+drm;

    for k =1:n2
        ux = atan(1/(xm2(k+1)*sqrt(1-(1/xm2(k+1))^2)));
        vx = Supersonic('nu','M',xm2(k+1),'g',gamma)*(pi/180.0);

        rxre(k+1) = c1*(1+0.5*gm1*xm2(k+1)*xm2(k+1));
        rxre(k+1) = rxre(k+1)^(gp1*0.5/gm1);
        rxre(k+1) = sqrt(abs(1-rxre(k+1)*sin(ve-vx+ux)/xp)); % abs is temporary

        xxre(k+1) = (1-rxre(k+1))*cos(ve-vx+ux)/sin(ve-vx+ux);

        xm2(k+2) = xm2(k+1) + drm;
        P2(k+2) = Supersonic('Po/P','M',xm2(k+2),'g',gamma);
    end
    P'
    P2'
%     plot(xx1re,rx1re,xx2re,rx2re,xxre,rxre)
    
    % my addition
    % cubic curve
%     cc(1,1) = xx1re(length(xx1re));
%     cc(1,2) = rx1re(length(rx1re));
%     cc(4,1) = xxre(1);
%     cc(4,2) = rxre(1);
%     dx1 = cc(4,1)-cc(1,1);
%     len = sqrt(((cc(4,1)-cc(1,1))^2)+((cc(4,2)-cc(1,2))^2));
%     a = rx1re(length(rx1re));
%     b = rx1re(length(rx1re)-1);
%     c = xx1re(length(xx1re));
%     d = xx1re(length(xx1re)-1);
%     ang1 = abs(atan((a-b)/(c-d)));
%     ang2 = abs(atan((rxre(2) - rxre(1))/(xxre(2) - xxre(1))));
%     cc(2,1) = cc(1,1) + ((len/3)*cos(ang1));
%     cc(2,2) = cc(1,2) - ((len/3)*sin(ang1));
%     cc(3,1) = cc(4,1) - ((len/3)*cos(ang2));
%     cc(3,2) = cc(4,2) + ((len/3)*sin(ang2));
% 
%     np = 25; % number of parameterization points
%     
%     for i = 1:np+1
%     	u(i) = (i-1)/np;
%     end
% 
%     n = length(cc); % degree of curve
%     b1 = decasteljau(u,n,np,cc);
%     
    xc = xx2re';
    yc = rx2re';
    xr = xx1re';
    yr = rx1re';

    for i = 1:length(xxre)
        xr(length(xx1re) + i) = xxre(i);
        yr(length(xx1re) + i) = rxre(i);
    end

%     
%     L1 = length(xr);
%     for i = 1:length(b1)
%         xr(L1+i) = b1(i,1);
%         yr(L1+i) = b1(i,2);
%     end
%     
%     L2 = length(xr);
%     for i = 1:length(xxre)
%         xr(L2+i) = xxre(i);
%         yr(L2+i) = rxre(i);
%     end
    
    scale = tlen/sqrt(((yc(1)-yr(1))^2)+((xc(1)-xr(1))^2));
    for i = 1:length(xc)
        xc(i) = xc(i)*scale;
        yc(i) = yc(i)*scale;
    end

    for i = 1:length(xr)
        xr(i) = xr(i)*scale;
        yr(i) = yr(i)*scale;
    end

    xdist = xr(1);
    ydist = yc(1);
    for i = 1:length(xc)
        xc(i) = xc(i) - xdist;
        yc(i) = yc(i) - ydist;
    end
    for i = 1:length(xr)
        xr(i) = xr(i) - xdist;
        yr(i) = yr(i) - ydist;
    end
end

% checks
function [] = check_linear_external(x, y)
    throat_len = sqrt((y(1)^2)+(x(1)^2))
    area_ratio = abs(y(length(y)))
end

function [] = check_linear_internal(xc, yc, xr, yr)
    throat_len = sqrt(((yc(1)-yr(1))^2)+((xc(1)-xr(1))^2))
%     throat_length = sqrt(((cowl(1,1)-ramp(1,1))^2)+((cowl(1,2)-ramp(1,2))^2))
%     throat_angle = atan((cowl(1,2)-ramp(1,2))/(cowl(1,1)-ramp(1,1)))*(180.0/pi)

    re = abs(yr(length(yr))) - abs(yc(length(yc)));
    area_ratio = re/throat_len
end

function [] = check_axisymmetric_external(x,y)
    re = abs(y(length(y)));
    r = abs(re - y(1));
    alpha = (pi/2)-atan((abs(x(1)))/(abs(y(1))));
    throat_len = abs((re-r)/sin(alpha))
    throat_len2 = sqrt((y(1)^2)+(x(1)^2))
    At = abs((pi/sin(alpha))*((re^2)-(r^2)));
    Ae = pi*(re^2);
    area_ratio = Ae/At
end

function [] = check_axisymmetric_internal(xc, yc, xr, yr)
    re = abs(yr(length(yr))) - abs(yc(length(yc)));
    Ae = pi*(re^2);
    
    re = abs(yr(length(yr)));
    r = re - abs(yr(1));
    alpha = (pi/2)-atan((abs(xr(1)))/(abs(yr(1))));
    At = abs((pi/sin(alpha))*((re^2)-(r^2)));

    throat_len = abs(re-r)/sin(alpha)
    throat_len2 = sqrt(((yc(1)-yr(1))^2)+((xc(1)-xr(1))^2))
    area_ratio = Ae/At
end
















