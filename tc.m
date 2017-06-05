syms r t

x = @(r,t) r*cos(t);
y = @(r,t) r*sin(t);

dxdr = diff(x,r);
dxdt = diff(x,t);
dydr = diff(y,r);
dydt = diff(y,t);
J = [dxdr dxdt; dydr dydt]

dxdrdr = diff(dxdr,r);
dxdrdt = diff(dxdr,t);
dxdtdr = diff(dxdt,r);
dxdtdt = diff(dxdt,t);
Hx = [dxdrdr dxdrdt; dxdtdr dxdtdt]

dydrdr = diff(dydr,r);
dydrdt = diff(dydr,t);
dydtdr = diff(dydt,r);
dydtdt = diff(dydt,t);
Hy = [dydrdr dydrdt; dydtdr dydtdt]

Jp = [cos(t) sin(t); -sin(t)/r cos(t)/r]


syms x y

r = @(x,y) sqrt((x^2)+(y^2));
t = @(x,y) atan(y/x);

drdx = simplify(diff(r,x));
drdy = simplify(diff(r,y));
dtdx = simplify(diff(t,x));
dtdy = simplify(diff(t,y));
Jp = [drdx drdy; dtdx dtdy]

drdxdx = simplify(diff(drdx,x));
drdxdy = simplify(diff(drdx,y));
drdydx = simplify(diff(drdy,x));
drdydy = simplify(diff(drdy,y));
Hr = [drdxdx drdxdy; drdydx drdydy]

dtdxdx = simplify(diff(dtdx,x));
dtdxdy = simplify(diff(dtdx,y));
dtdydx = simplify(diff(dtdy,x));
dtdydy = simplify(diff(dtdy,y));
Ht = [dtdxdx dtdxdy; dtdydx dtdydy]










