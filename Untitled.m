clear all
clc

% curve fit
x = [40.5585, 63.5331, 87.5466, ...
     44.2217, 68.5106, 93.3247, ...
     47.2535, 71.9528, 97.3989]'
 
y = [9.9302, 14.8313, 19.7762, ...
     9.8586, 14.8037, 19.7112, ...
     9.9147, 14.8001, 19.7110]'

k = 1;
for r = 1:4:9
    for A = 10:5:20    
        C(k,1) = A*A*r*r;
        C(k,2) = A*A*r;
        C(k,3) = A*A;
        C(k,4) = A*r*r;
        C(k,5) = A*r;
        C(k,6) = A;
        C(k,7) = r*r;
        C(k,8) = r;
        C(k,9) = 1.0;
        
        k = k + 1;
    end
end

Cx = C\x
Cy = C\y

myfuncx = @(x,y) ((Cx(1)*x*x*y*y)+(Cx(2)*x*x*y)+(Cx(3)*x*x)+ ...
                 (Cx(4)*x*y*y)+(Cx(5)*x*y)+(Cx(6)*x)+ ...
                 (Cx(7)*y*y)+(Cx(8)*y)+Cx(9));
             
myfuncy = @(x,y) ((Cy(1)*x*x*y*y)+(Cy(2)*x*x*y)+(Cy(3)*x*x)+ ...
                 (Cy(4)*x*y*y)+(Cy(5)*x*y)+(Cy(6)*x)+ ...
                 (Cy(7)*y*y)+(Cy(8)*y)+Cy(9));

X = myfuncx(12.375, 2.48)
Y = myfuncy(12.375, 2.48)             

% newtons method
myx = 53.086285;
myy = 12.7;

dfxdA = @(A,r) (2.0*Cx(1)*A*r*r)+(2.0*Cx(2)*A*r)+(2.0*Cx(3)*A)+(Cx(4)*r*r)+(Cx(5)*r)+Cx(6);
dfxdr = @(A,r) (2.0*Cx(1)*A*A*r)+(Cx(2)*A*A)+(2.0*Cx(4)*A*r)+(Cx(5)*A)+(2.0*Cx(7)*r)+Cx(8);

dfydA = @(A,r) (2.0*Cy(1)*A*r*r)+(2.0*Cy(2)*A*r)+(2.0*Cy(3)*A)+(Cy(4)*r*r)+(Cy(5)*r)+Cy(6);
dfydr = @(A,r) (2.0*Cy(1)*A*A*r)+(Cy(2)*A*A)+(2.0*Cy(4)*A*r)+(Cy(5)*A)+(2.0*Cy(7)*r)+Cy(8);

x = [15;2]
for i = 1:100
F = [myfuncx(x(1),x(2)) - myx; ...
    myfuncy(x(1),x(2)) - myy]

J = [dfxdA(x(1),x(2)), ...
     dfxdr(x(1),x(2)); ...
     dfydA(x(1),x(2)), ...
     dfydr(x(1),x(2))]
 inv(J)
 y = -J\F
 x = x + y
end