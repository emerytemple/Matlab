function [x, y] = Untitled(xp, tlen, gama, nr) % axi-symmetric external contour
% xp - expansion ratio
% tlen - length of throat
% gama - ratio of specific heats
% nr - number of discretization points
%
% throat is from (0,0) to (x(1),y(1))
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
    
    plot(x,y)
end

function [] = Untitled2() % axi-symmetric internal-external contour
    clear all
    clc
    format long
    
    n1 = 50
    n2 = 50

    xp = 6.5
    rrre = 0.3
    pht = pi/2.0

    gama = 1.4
    gp1 = gama+1.0;
    gm1 = gama-1.0;

%     peipc = 0.0994474;
%     rmei = peipc^((1.0-gama)/gama);
%     rmei = (2.0/gm1)*(rmei-1.0);
%     rmei = sqrt(rmei);
%     xpei = Supersonic('A/A*','M',rmei,'g',gama)
    xpei = 1.7
    
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