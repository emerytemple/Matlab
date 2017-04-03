function [] = throat_test()
%     M = 3.2;
%     TL = 0.01;
% 
%     y1 = -Supersonic('A/A*','M',1.0,'g',1.4)*TL*sind(Supersonic('mu','M',1.0)+Supersonic('nu','M',M,'g',1.4));
%     y2 = -M*Supersonic('A/A*','M',M,'g',1.4)*TL*sind(Supersonic('mu','M',M));
%     yt = y1-y2-4;
%     
%     fff = @(M,TL) ((-Supersonic('A/A*','M',1.0,'g',1.4)*TL*sind(Supersonic('mu','M',1.0)+Supersonic('nu','M',M,'g',1.4)))-(-M*Supersonic('A/A*','M',M,'g',1.4)*TL*sind(Supersonic('mu','M',M)))-4);
%     
%     fff(3.2,1)

    hold on
    j = 1;
    for i = 1.4:0.25:4.4
        fff = @(M,TL) ((-Supersonic('A/A*','M',1.0,'g',1.4)*TL*sind(Supersonic('mu','M',1.0)+Supersonic('nu','M',M,'g',1.4)))-(-M*Supersonic('A/A*','M',M,'g',1.4)*TL*sind(Supersonic('mu','M',M)))-8);
        ansss = fzero(@(x) fff(i,x), [0.01,100.0]);
    
        [xx,yy] = LEA(i,ansss,1.4,20);
        th1x = [0,-xx(1)];
        th1y = [0,-yy(1)];
    
        bbb(j,1) = -xx(1);
        bbb(j,2) = -yy(1);
        j = j + 1;
        
        plot(th1x,th1y,xx-xx(1),yy-yy(1))    
    end
    plot(bbb(:,1),bbb(:,2))
    hold off
    bbb
    yy(length(yy))-yy(1)
    
%     hold on
%     for i = 2.2:1:4.2
%         [xx,yy] = LEA(i, 3, 1.4, 20);
% 
%         th1x = [0,xx(1)];
%         th1y = [0,yy(1)];
% 
%         plot(th1x,th1y,xx,yy) 
%     end
%     hold off

end

function [x,y] = LEA(mach_exit, throat_len, gamma, nr) % linear external (approximate) contour
% Angelino, G. - Approximate Method for Plug Nozzle Design (AIAA)

    x = zeros(nr,1);
    y = zeros(nr,1);
    
    dm = (mach_exit-1.0)/(nr-1.0);
    theta = Supersonic('nu','M',mach_exit,'g',gamma); % ve

    for i = 1:nr
        M = 1.0+((i-1)*dm);
        mu = Supersonic('mu','M',M);
        nu = Supersonic('nu','M',M,'g',gamma);
        alpha = mu+theta-nu;
        xsi = M*Supersonic('A/A*','M',M,'g',gamma);

        % redimensionalize
        L = xsi*throat_len;
        x(i) = L*cosd(alpha);
        y(i) = -L*sind(alpha);
    end
end