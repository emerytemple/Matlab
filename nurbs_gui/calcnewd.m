function [dx,dy] = calcnewd(cp,op)
    p1x = cp(2,1);
    p1y = cp(2,2);
    p2x = cp(3,1);
    p2y = cp(3,2);

    m1 = (p2y-p1y)/(p2x-p1x);
    b1 = p1y-(m1*p1x);

    p3x = op(2,1);
    p3y = op(2,2);
    p4x = op(3,1);
    p4y = op(3,2);

    m2 = (p4y-p3y)/(p4x-p3x);
    b2 = p3y-(m2*p3x);

    dx = (b2-b1)/(m1-m2);
    dy = (m1*dx)+b1;
end
