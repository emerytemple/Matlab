function [] = Untitled()
    clear all
    clc
    format long
    
    n1 = 50
    n2 = 50

    xp = 12.0
    rrre = 0.3
    pht = 1.570796
    
    gama = 1.2
    gp1 = gama+1.0;
    gm1 = gama-1.0;

    peipc = 0.0994474;
    rmei = peipc^((1.0-gama)/gama);
    rmei = (2.0/gm1)*(rmei-1.0);
    rmei = sqrt(rmei);
    xpei = Supersonic('A/A*','M',rmei,'g',gama)
    
    rm = Supersonic('Mar','x',xp,'g',gama);
    ve = Supersonic('nu','M',rm,'g',gama)*(pi/180.0);
    rmei = Supersonic('Mar','x',xpei,'g',gama);
    vei = Supersonic('nu','M',rmei,'g',gama)*(pi/180.0);
    y = 1.0/rmei;
    uei = atan(y/sqrt(1.0-y*y));
    thei = vei - ve;
    phei = thei + uei;
    ToT = Supersonic('To/T','M',rmei,'g',gama);
    a1 = (2.0/gp1)*ToT;
    a1 = a1^(gp1/(2.0*gm1));
    b1 = sin(phei);
    rpre = sqrt(1.0-a1*b1/xp);
    xpre = (rpre-1.0)*cos(phei)/sin(phei);
    drm = (rmei-1.0)/n1;
    
    % internal section
    k = 0;
    xm = 1.0;
    while((k-n1) <= 0.0)
        vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
        bx = pht-1.570796-vx+abs(thei);
        x1re = 2.0*rrre*sin(0.5*bx);
        psi = 3.1416-pht+vx-0.5*(3.1416-bx);
        rx1re = rpre+x1re*sin(psi);
        xx1re = xpre-x1re*cos(psi);
        if k == 0
            rx2re = sqrt(rx1re*rx1re+sin(pht)/xp);
            xx2re = xx1re+(rx2re-rx1re)*cos(pht)/sin(pht);
        else
            ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
            phx = 2.0*vei-ve-vx+ux;
            ToT = Supersonic('To/T','M',xm,'g',gama);
            a2 = (2.0/gp1)*ToT;
            b2 = 0.5*gp1/gm1;

            rx2re = sqrt(rx1re*rx1re+(a2^b2)*sin(phx)/xp);
            xx2re = xx1re+(rx2re-rx1re)*cos(phx)/sin(phx);
        end
        
        ansss(k+1,1) = xm;
        ansss(k+1,2) = rx1re;
        ansss(k+1,3) = xx1re;
        ansss(k+1,4) = rx2re;
        ansss(k+1,5) = xx2re;
        
        k = k+1;
        xm = xm+drm;
    end
    xm = xm-drm;
    ansss
    
    % external section
    ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
    vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
    ToT = Supersonic('To/T','M',xm,'g',gama);
    rxre = (2.0/gp1)*ToT;
    rxre = rxre^(gp1/(2.0*gm1));
    rxre = 1.0-rxre*sin(ve-vx+ux)/xp;
    rxre = sqrt(rxre);
    xxre = (1.0-rxre)*cos(ve-vx+ux)/sin(ve-vx+ux);
    
	ret(1,1) = xm;
    ret(1,2) = rxre;
    ret(1,3) = xxre; 
    
    k1 = 1;
    drm = (rm-xm)/n2;
    xm = xm+drm;
    while((k1-n2) <= 0)
        ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
        vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
        ToT = Supersonic('To/T','M',xm,'g',gama);
        rxre = (2.0/gp1)*ToT;
        rxre = rxre^(gp1*0.5/gm1);
        rxre = sqrt(1.0-rxre*sin(ve-vx+ux)/xp);
        xxre = (1.0-rxre)*cos(ve-vx+ux)/sin(ve-vx+ux);

        ret(k1+1,1) = real(xm);
        ret(k1+1,2) = real(rxre);
        ret(k1+1,3) = real(xxre);

        xm = xm+drm;
        k1 = k1+1;
    end
    ret
    
    plot(ansss(:,3),ansss(:,2),ansss(:,5),ansss(:,4),ret(:,3),ret(:,2))
end

function [] = mycontour(xp, tlen, tang, gamma, n1, n2, rrre)
% xp - end of external area ratio
% tang - angle (in degrees) the throat differs from vertical

    xpe = 75.13;
    tang = 0.0;
    re = 1;
    n2 = 30;

    Mt = 1.0
    vt = 0.0
    mut = 90.0
    
    Me = Supersonic('Mar','x',xpe,'g',gamma)
    ve = Supersonic('nu','M',Me,'g',gamma)
    mue = Supersonic('mu','M',Me)
    
    vei = 0.5*(ve-tang)
    Mei = Supersonic('Mnu','v',vei,'g',gamma)
    muei = Supersonic('mu','M',Mei)
    
    Ae = pi*re*re
    At = Ae/xpe
    
    th3 = ve-vei
    th2 = vei-vt
    
    % external section
    dth = th3/(n2-1);
    for i = 1:n2
        th = th3 - ((i-1)*dth);
        if th < 0
            th = 0;
        end
        v = vei + th;
        M = Supersonic('Mnu','v',v,'g',gamma);
        mu = Supersonic('mu','M',M);
        xp = Supersonic('A/A*','M',M,'g',gamma);
        A = xp*At;
        r = sqrt(A/pi);
        P(i,1) = r/tand(mu);
        P(i,2) = -r;
    end
    P
    plot(P(:,1),P(:,2))
end

function [xc, yc, xr, yr] = axisymmetric_internal(xp, tlen, tang, gamma, n1, n2, rrre) % 
% axi-symmetric internal-external contour of lee & thompson
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

function [varargout] = Supersonic(string,varargin)

% set inputs
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'b' % angle of oblique shock
            beta1 = varargin{i+1};
            beta1 = beta1*(pi/180);
        case 'g' % ratio of specific heats
            gama = varargin{i+1};
            gp1 = gama+1.0;
            gm1 = gama-1.0;
        case 'M' % mach number
            M = varargin{i+1};
        case 'p' % specific heat
            cp = varargin{i+1};
%         case 'P' % pressure
%             P = varargin{i+1};
        case 'PR' % pressure ratio (Po/P)
            pr = varargin{i+1};
        case 'r' % density
            rho = varargin{i+1};
        case 'R' % ideal gas law constant
            R = varargin{i+1};
        case 't' % new flow angle for oblique shock
            theta = varargin{i+1};
            theta = theta*(pi/180);
%         case 'T' % temperature
%             T = varargin{i+1};
        case 'v' % prandtl-meyer function
            nu = varargin{i+1};
        case 'x' % area ratio
            xp = varargin{i+1};
        otherwise
            error = 1
    end
end

% calculate value
switch string

    % 1-D
    case 'P*/Po'
        varargout{1} = (2.0/gp1)^(gama/gm1);
    case 'Po/P'
        ToT = 1.0+((gm1/2.0)*(M^2));
        varargout{1} = ToT^(gama/gm1);    
    case 'r*/ro'
        varargout{1} = (2.0/gp1)^(1.0/gm1);
    case 'ro/r'
        ToT = 1.0+((gm1/2.0)*(M^2));
        varargout{1} = ToT^(1.0/gm1);
    case 'T*/To'
        varargout{1} = 2.0/gp1;
    case 'To/T'
        varargout{1} = 1.0+((gm1/2.0)*(M^2));

    case 'Mpr'
        fmf = @(m) (pr-Supersonic('Po/P','M',m,'g',gama));
        varargout{1} = fzero(@(x) fmf(x), [0.01,50.0]);
        
    % 1-D w/ heat addition
    case 'P/P*h'
        varargout{1} = gp1/(1.0+(gama*(M^2)));
    case 'T/T*h'
        PPs = gp1/(1.0+(gama*(M^2)));
        varargout{1} = (M^2)*(PPs^2);
    case 'r/r*h'
        PPs = gp1/(1.0+(gama*(M^2)));
        varargout{1} = (1/(M^2))*(1/PPs);
    case 'Po/Po*h'
        PPs = gp1/(1.0+(gama*(M^2)));
        varargout{1} = PPs*(((2.0+(gm1*(M^2)))/gp1)^(gama/gm1));
    case 'To/To*h'
        varargout{1} = ((gp1*(M^2))/((1.0+(gama*(M^2)))^2))*(2.0+(gm1*(M^2)));
        
    % 1-D w/ friction
    case 'T/T*f'
        varargout{1} = gp1/(2.0+(gm1*(M^2)));
    case 'P/P*f'
        TTs = gp1/(2.0+(gm1*(M^2)));
        varargout{1} = sqrt(TTs)/M;
    case 'r/r*f'
        TTs = gp1/(2.0+(gm1*(M^2)));
        varargout{1} = sqrt(1/TTs)/M;
    case 'Po/Po*f'
        TTs = gp1/(2.0+(gm1*(M^2)));
        varargout{1} = ((1/TTs)^(gp1/(2.0*gm1)))/M;
    case '4fL*/D'
        TTs = gp1/(2.0+(gm1*(M^2)));
        varargout{1} = ((1-(M^2))/(gama*(M^2)))+((gp1/(2.0*gama))*log(TTs*(M^2)));
        
    % quasi 1-D
	case 'A/A*'
        ToT = 1.0+((gm1/2.0)*(M^2));
        varargout{1} = sqrt((1/(M^2))*(((2.0/gp1)*ToT)^(gp1/gm1)));
    case 'Msub'
        % xp <= 1.0
        gp1o2gm1 = gp1/(2*gm1);
        Togp1 = 2/gp1;
        gm1o2 = gm1/2;

        fmf = @(M,xp,gama) ((M*xp)-((Togp1*(1+(gm1o2*(M^2))))^gp1o2gm1));
        varargout{1} = fzero(@(x) fmf(x,xp,gama), [0.0,1.0]);
    case 'Msup'
        % xp >= 1.0
        gp1o2gm1 = gp1/(2*gm1);
        Togp1 = 2/gp1;
        gm1o2 = gm1/2;

        fmf = @(M,xp,gama) ((M*xp)-((Togp1*(1+(gm1o2*(M^2))))^gp1o2gm1));        
        varargout{1} = fzero(@(x) fmf(x,xp,gama), [1.0,50.0]);
    case 'Mar'
        if xp == 1.0
            varargout{1} = 1.0;
        else
            fmf = @(m) (xp-Supersonic('A/A*','M',m,'g',gama));
            varargout{1} = fzero(@(x) fmf(x), [1.0,50.0]);
        end
        
    % normal shock
    case 'M2'
        % M >= 1.0
        ToT = 1.0+((gm1/2.0)*(M^2));
        varargout{1} = sqrt(ToT/((gama*(M^2))-(gm1/2.0)));
    case 'Po2/Po1'
        PoP = 1.0+(((2.0*gama)/gp1)*((M^2)-1.0));
        ror = (gp1*(M^2))/(2.0+(gm1*(M^2)));
        ToT = PoP/ror;
        del_s = (cp*log(ToT))-(R*log(PoP));
        varargout{1} = exp(-del_s/R);
    case 'Po2/P1'
        PoP = 1.0+(((2.0*gama)/gp1)*((M^2)-1.0));
        ror = (gp1*(M^2))/(2.0+(gm1*(M^2)));
        ToT = PoP/ror;
        del_s = (cp*log(ToT))-(R*log(PoP));
        Po2Po1 = exp(-del_s/R);
        TT = 1.0+((gm1/2.0)*(M^2));
        Po1P1 = (TT)^(gama/(gama-1.0));
        varargout{1} = Po2Po1*Po1P1; 
    case 'P2/P1'
        varargout{1} = 1.0+(((2.0*gama)/gp1)*((M^2)-1.0));
    case 'r2/r1'
        varargout{1} = (gp1*(M^2))/(2.0+(gm1*(M^2)));
    case 'T2/T1'
        PoP = 1.0+(((2.0*gama)/gp1)*((M^2)-1.0));
        ror = (gp1*(M^2))/(2.0+(gm1*(M^2)));
        varargout{1} = PoP/ror;
        
    % oblique shock wave    
    case 'theta'
        tcb = 2.0*cot(beta1);
        top = ((M*sin(beta1))^2)-1.0;
        bot = ((M^2)*(gama+cos(2.0*beta1)))+2.0;
        theta = atan(tcb*(top/bot));
        varargout{1} = theta*(180/pi);
    case 'betas'
        m2 = ((M^2)-1.0)^2;
        m3 = ((M^2)-1.0)^3;
        gm2 = (gm1/2.0)*(M^2);
        gp2 = (gp1/2.0)*(M^2);
        gp4 = (gp1/4.0)*(M^4);
        t2 = tan(theta)^2;
        lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2));
        xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/(lam^3);
        delta = 1.0;
        cx = acos(xi);
        part = ((4.0*pi*delta)+cx)/3.0;
        cp = cos(part);
        tb = ((M^2)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta));
        beta1 = atan(tb);
        varargout{1} = beta1*(180/pi);
    case 'betaw'
        m2 = ((M^2)-1.0)^2;
        m3 = ((M^2)-1.0)^3;
        gm2 = (gm1/2.0)*(M^2);
        gp2 = (gp1/2.0)*(M^2);
        gp4 = (gp1/4.0)*(M^4);
        t2 = tan(theta)^2;
        lam = sqrt(m2-(3.0*(1.0+gm2)*(1.0+gp2)*t2));
        xi = (m3-(9.0*(1.0+gm2)*(1.0+gm2+gp4)*t2))/(lam^3);
        delta = 0.0;
        cx = acos(xi);
        part = ((4.0*pi*delta)+cx)/3.0;
        cp = cos(part);
        tb = ((M^2)-1.0+(2.0*lam*cp))/(3.0*(1.0+gm2)*tan(theta));
        beta1 = atan(tb);
        varargout{1} = beta1*(180/pi);
    case 'Mob'
        ftbm = @(t,b,m) (tan(t) - ...
                        (2.0*cot(b)* ...
                        ((((m*sin(b))^2)-1.0)/ ...
                        (((m^2)*(gama+cos(2.0*b)))+2.0))));
        varargout{1} = fzero(@(x) ftbm(theta,beta1,x),[1.0,50.0]);
        
    % prandtl-meyer expansion waves    
	case 'mu'
        varargout{1} = asind(1/M);
    case 'nu'
        m2 = (M^2)-1.0;
        varargout{1} = (sqrt(gp1/gm1)*atand(sqrt((gm1/gp1)*m2)))-atand(sqrt(m2));
    case 'Mnu'
        gp1ogm1 = gp1/gm1;
        gm1ogp1 = gm1/gp1;

        nu = nu*(pi/180);

        fpmf = @(M) ((sqrt(gp1ogm1)*atan(sqrt(gm1ogp1*((M^2)-1))))-atan(sqrt((M^2)-1))-nu);

        varargout{1} = fzero(fpmf, [1,50]);
    otherwise
        error = 1
end

end





























%     n1 = 50
%     n2 = 50
%     
%     gas = 1.0
%     peipc = 0.0994474
%     xp = 12.0
%     rrre = 0.3
%     rm = 4.1
%     pht = 1.570796
%     
%     r = 1718.0
%     te = 3000.0
%     papc = 0.01
%     g = 32.17405
%     
%     gama = 1.2
%     
%     rm = Supersonic('Mar','x',xp,'g',gama);
%     ve = Supersonic('nu','M',rm,'g',gama)*(pi/180.0);
% 	rmei = peipc^((1.0-gama)/gama);
%     rmei = (2.0/(gama-1.0))*(rmei-1.0);
%     rmei = sqrt(rmei);
%     vei = Supersonic('nu','M',rmei,'g',gama)*(pi/180.0);
%     y = 1.0/rmei;
%     uei = atan(y/sqrt(1.0-y*y));
%     thei = vei - ve;
%     phei = thei + uei;
%     a1 = (2.0/(gama+1.0))*(1.0+0.5*(gama-1.0)*rmei*rmei);
%     a1 = a1^((gama+1.0)/(2.0*(gama-1.0)));
%     b1 = sin(phei);
%     rpre = sqrt(1.0-a1*b1/xp);
%     xpre = (rpre-1.0)*cos(phei)/sin(phei);
%     drm = (rmei-1.0)/n1;
%     
%     % internal section
%     k = 0;
%     xm = 1.0;
%     while((k-n1) <= 0.0)
%         vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
%         bx = pht-1.570796-vx+abs(thei);
%         x1re = 2.0*rrre*sin(0.5*bx);
%         psi = 3.1416-pht+vx-0.5*(3.1416-bx);
%         rx1re = rpre+x1re*sin(psi);
%         xx1re = xpre-x1re*cos(psi);
%         if k == 0
%             rx2re = sqrt(rx1re*rx1re+sin(pht)/xp);
%             xx2re = xx1re+(rx2re-rx1re)*cos(pht)/sin(pht);
%         else
%             ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
%             phx = 2.0*vei-ve-vx+ux;
%             a2 = (2.0/(gama+1.0))*(1.0+0.5*(gama-1.0)*xm*xm);
%             b2 = 0.5*(gama+1.0)/(gama-1.0);
% 
%             rx2re = sqrt(rx1re*rx1re+(a2^b2)*sin(phx)/xp);
%             xx2re = xx1re+(rx2re-rx1re)*cos(phx)/sin(phx);
%         end
%         pxpc = (1.0+0.5*(gama-1.0)*xm*xm)^(-gama/(gama-1.0));
%         
%         ansss(k+1,1) = xm;
%         ansss(k+1,2) = rx1re;
%         ansss(k+1,3) = xx1re;
%         ansss(k+1,4) = rx2re;
%         ansss(k+1,5) = xx2re;
%         ansss(k+1,6) = pxpc;
%         
%         k = k+1;
%         xm = xm+drm;
%     end
%     xm = xm-drm;
%     ansss
%     
%     % external section
%     ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
%     vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
%     rxre = (2.0/(gama+1.0))*(1.0+0.5*(gama-1.0)*xm*xm);
%     rxre = rxre^((gama+1.0)/(2.0*(gama-1.0)));
%     rxre = 1.0-rxre*sin(ve-vx+ux)/xp;
%     rxre = sqrt(rxre);
%     xxre = (1.0-rxre)*cos(ve-vx+ux)/sin(ve-vx+ux);
%     c1 = (2.0/(gama+1.0))^(gama/(gama-1.0));
%     c2 = sqrt((0.5*(gama+1.0)*xm*xm)/(1.0+0.5*(gama-1.0)*xm*xm));
%     sumcg = gama*c1*c2*cos(thei)+xp*pxpc*(1.0-rxre*rxre);
%     vt = te*(1.0+0.5*(gama-1.0)*rm*rm);
%     vt = gama*g*r*vt/(0.5*(gama+1.0));
%     vt = sqrt(vt);
%     vq = te*(1.0+0.5*(gama-1.0)*rm*rm);
%     vc = 1.0+0.5*(gama-1.0)*xm*xm;
%     vq = gama*r*g*vq/vc;
%     vq = sqrt(vq);
%     a = 1.0+0.5*(gama-1.0)*xm*xm;
%     b = -gama/(gama-1.0);
%     a = a^b;
%     c = 1.0-papc*a*sin(phei)/(xm*xm);
%     d = cos(thei)+c/gama;
%     sumim = vq*d/g;
%     co = cos(thei)+1.0/gama;
%     sumva = vq*co/g;
%     
% 	ret(1,1) = xm;
%     ret(1,2) = rxre;
%     ret(1,3) = xxre;
%     ret(1,4) = pxpc;
%     ret(1,5) = sumcg;
%     ret(1,6) = sumim;
%     ret(1,7) = sumva;    
%     
%     k1 = 1;
%     drm = (rm-xm)/n2;
%     xm = xm+drm;
%     pro = pxpc;
%     rxo = rxre;
%     while((k1-n2) <= 0)
%         ux = atan(1.0/(xm*sqrt(1.0-(1.0/xm)^2)));
%         vx = Supersonic('nu','M',xm,'g',gama)*(pi/180.0);
%         rxre = (2.0/(gama+1.0))*(1.0+0.5*(gama-1.0)*xm*xm); %
%         rxre = rxre^((gama+1.0)*0.5/(gama-1.0));
%         rxre = sqrt(1.0-rxre*sin(ve-vx+ux)/xp);
%         xxre = (1.0-rxre)*cos(ve-vx+ux)/sin(ve-vx+ux);
%         pxpc = (1.0+0.5*(gama-1.0)*xm*xm)^(-gama/(gama-1.0));
%         sumcg = sumcg+0.5*xp*(pro+pxpc)*(rxo*rxo-rxre*rxre);
%         a = gama/(gama-1.0);
%         a = (0.5*(gama+1.0))^a;
%         a = a*vt/(gama*g);
%         b = 0.5*a*xp;
%         sumim = sumim+b*(pro+pxpc-2.0*papc)*(rxo*rxo-rxre*rxre);
%         sumva = sumva+b*(pro+pxpc)*(rxo*rxo-rxre*rxre);
% 
%         ret(k1+1,1) = real(xm);
%         ret(k1+1,2) = real(rxre);
%         ret(k1+1,3) = real(xxre);
%         ret(k1+1,4) = real(pxpc);
%         ret(k1+1,5) = real(sumcg);
%         ret(k1+1,6) = real(sumim);
%         ret(k1+1,7) = real(sumva);
% 
%         xm = xm+drm;
%         pro = pxpc;
%         rxo = rxre;
%         k1 = k1+1;
%     end
%     ret
%     
%     plot(ansss(:,3),ansss(:,2),ansss(:,5),ansss(:,4),ret(:,3),ret(:,2))
% end
% 
% function [] = mycontour(xp, tlen, tang, gamma, n1, n2, rrre)
% % xp - end of external area ratio
% % tang - angle (in degrees) the throat differs from vertical
% 
%     xpe = 75.13;
%     tang = 0.0;
%     re = 1;
%     n2 = 30;
% 
%     Mt = 1.0
%     vt = 0.0
%     mut = 90.0
%     
%     Me = Supersonic('Mar','x',xpe,'g',gamma)
%     ve = Supersonic('nu','M',Me,'g',gamma)
%     mue = Supersonic('mu','M',Me)
%     
%     vei = 0.5*(ve-tang)
%     Mei = Supersonic('Mnu','v',vei,'g',gamma)
%     muei = Supersonic('mu','M',Mei)
%     
%     Ae = pi*re*re
%     At = Ae/xpe
%     
%     th3 = ve-vei
%     th2 = vei-vt
%     
%     % external section
%     dth = th3/(n2-1);
%     for i = 1:n2
%         th = th3 - ((i-1)*dth);
%         if th < 0
%             th = 0;
%         end
%         v = vei + th;
%         M = Supersonic('Mnu','v',v,'g',gamma);
%         mu = Supersonic('mu','M',M);
%         xp = Supersonic('A/A*','M',M,'g',gamma);
%         A = xp*At;
%         r = sqrt(A/pi);
%         P(i,1) = r/tand(mu);
%         P(i,2) = -r;
%     end
%     P
%     plot(P(:,1),P(:,2))
