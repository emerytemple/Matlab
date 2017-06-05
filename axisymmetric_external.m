function [x,y] = axisymmetric_external(xpi, xp, re, gamma, nr) % axi-symmetric external contour
    x = zeros(nr,1);
    y = zeros(nr,1);

    mach_exit = Supersonic('Mar','x',xp,'g',gamma);
    dm = (mach_exit-1.0)/(nr-1.0);

    theta = Supersonic('nu','M',mach_exit,'g',gamma) % ve
    theta_rad = theta*(pi/180.0);

%     re = xp*tlen*(1+sqrt(1-(cos(theta_rad)/xp)));
    
    % calculate contour
    for i = 1:nr
        M = 1.0+((i-1)*dm);
        mu = Supersonic('mu','M',M);
        nu = Supersonic('nu','M',M,'g',gamma);
        alpha = (mu+theta-nu)*(pi/180.0);
        eps = Supersonic('A/A*','M',M,'g',gamma);
        xsi = (1.0-sqrt(1.0-((eps/xp)*M*sin(alpha))))/sin(alpha);
        
        % redimensionalize
        L = real(xsi)*re;
        x(i) = L*cos(alpha);
        y(i) = -L*sin(alpha);
        
%         ansss(i,1) = i;
%         ansss(i,2) = x(i);
%         ansss(i,3) = y(i);
    end
    
%     ansss
    
%     check_axisymmetric_external(x,y)

    % find split point
    M = Supersonic('Mar','x',xpi,'g',gamma);
    mu = Supersonic('mu','M',M);
    nu = Supersonic('nu','M',M,'g',gamma);
    alpha = (mu+theta-nu)*(pi/180.0);
    eps = Supersonic('A/A*','M',M,'g',gamma);
    xsi = (1.0-sqrt(1.0-((eps/xp)*M*sin(alpha))))/sin(alpha);
    L = xsi*re;
    xc = L*cos(alpha)
    yc = -L*sin(alpha)
    
    % calculate initial point
    mach_inlet = Supersonic('Mar','x',xpi,'g',gamma)
    theta = Supersonic('nu','M',mach_exit,'g',gamma) - Supersonic('nu','M',mach_inlet,'g',gamma)
    theta_rad = theta*(pi/180.0)
    
    % alternate method
%     Ae = pi*re*re
%     A = (xpi/xp)*Ae
%     r = sqrt((re*re)-((A*cos(theta_rad))/pi))
%     el = r/cos(theta_rad)
%     L = re/cos(theta_rad)
%     len = L-el
%     xp1 = -len*sin(theta_rad)
%     yp1 = -len*cos(theta_rad)
%     theta_check = atan(yp1/xp1)
%     area_check = pi
    
    % find intersection of lines
    xp = (tan(theta_rad)*(yc+(tan(theta_rad)*xc)))/(sec(theta_rad)^2)
    yp = xp/tan(theta_rad)
    
    % split contour
    for i = 1:nr
        if xc < x(i)
            ind = i;
            break;
        end
    end

    % generate acutal contour
    final(1,1) = xp;
    final(1,2) = yp;
    
    final(2,1) = xc;
    final(2,2) = yc;
    
    for i = ind:nr
        final(i-ind+3,1) = x(i);
        final(i-ind+3,2) = y(i);
    end
    
    
    final
    
    plot(final(:,1),final(:,2),xp,yp,'o',xc,yc,'o',x,y)
    check_axisymmetric_external(final(:,1),final(:,2))
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