function [] = cfdtest()

clear all
clc

x = 0:0.05:10;
a = 1; % m/s
tf = 0.75; % sec

res2 = fdiff(x);

% plot(x, res1, x, res2)
% xlim([-0.1, 5])
% ylim([-0.1, 1.1])

end

function [y] = myf(x)
    if x < 0.0
        y = 1.0;
    else
        y = 0.0;
    end
end

function [res] = actual(x, a, t)
    for i = 1:length(x)
        res(i) = myf(x(i) - (a*t));
    end
end

function [res] = fdiff(x)
    cfl = 0.8;
    a = 1;
    dx = x(2) - x(1);
    dt = (cfl*dx)/a;
    
    t = -dt;
    
    u(:,1) = actual(x, a, 0.5);
    plot(x, actual(x, a, 0.5), x, u(:,1))
    xlim([-0.1, 10.1])
    ylim([-0.1, 1.1])
    F(1) = getframe;
    
    for n = 1:length(x)
        t = t + dt;
        u(1,n+1) = 1.0;
        u(2,n+1) = 1.0;
        for i = 3:length(x)-2
            % ftcs (unconditionally unstable) (2nd order)
            % u(i,n+1) = u(i,n) - (0.5*cfl*(u(i+1,n)-u(i-1,n)));
            
            % ftbs (conditionally stable) (1st order)
            % u(i,n+1) = u(i,n) - (cfl*(u(i,n)-u(i-1,n)));
            
            % lax-friedrichs
            % u(i,n+1) = (0.5*(u(i+1,n)+u(i-1,n))) - (0.5*cfl*(u(i+1,n)-u(i-1,n)));
            
            % lax-wendroff
            % u(i,n+1) = u(i,n) - (0.5*cfl*(u(i+1,n)-u(i-1,n))) + (0.5*cfl*cfl*(u(i+1,n)-(2.0*u(i,n))+u(i-1,n)));
            
            % lax-wendroff w/ superbee limiter
            if abs(u(i+1,n) - u(i,n)) < 0.00001
                ri = 0.0;
            else
                ri = (u(i,n) - u(i-1,n)) / (u(i+1,n) - u(i,n));
            end
            if abs(u(i,n) - u(i-1,n)) < 0.00001
                rim1 = 0.0;
            else
                rim1 = (u(i-1,n)-u(i-2,n)) / (u(i,n)-u(i-1,n));
            end
            p1 = cfl*(u(i,n) - u(i-1,n));
            p2 = 0.5*cfl*(1.0-cfl)*xsii(ri)*(u(i+1,n)-u(i,n));
            p3 = 0.5*cfl*(1.0-cfl)*xsii(rim1)*(u(i,n)-u(i-1,n));
            u(i,n+1) = u(i,n) - p1 - p2 + p3;
        end
        u(length(x),n+1) = 0.0;
        u(length(x),n+1) = 0.0;
        
        plot(x, actual(x, a, t+10*dt), x, u(:,n+1))
        xlim([-0.1, 10.1])
        ylim([-0.1, 1.1])
        F(n+1) = getframe;
    end
    res = u(:,length(x));
    movie(F)
end

function [res] = xsii(r)
    res = max([0.0, min(2.0*r, 1.0), min(r, 2.0)]);
end








