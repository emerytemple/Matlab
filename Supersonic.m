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