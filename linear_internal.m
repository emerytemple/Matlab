function [xc, yc, xr, yr] = linear_internal(tht, xpe, rrre, gamma, nei, ne)
%     gamma = 1.4;
%     xpe = 5.0;
%     tht = 0.0;
%     rt = 1.0;
%     rrre = 0.5;
%     nei = 48;
%     ne = 47;
    rt = 1.0; % non-dimensional contour

    me = Supersonic('Mar','x',xpe,'g',gamma);
    ve = Supersonic('nu','M',me,'g',gamma);
    vei = (ve-tht)/2.0;

    tht = tht*(pi/180.0);
    vei = vei*(pi/180.0);
    ve = ve*(pi/180.0);

    xc(1) = 0.0;
    yc(1) = 0.0;

    xr(1) = -rt*sin(tht);
    yr(1) = rt*cos(tht);

    thx = tht;
    thxc = tht;

    dv = (vei-0.0)/(nei-1.0);
	for i = 2:nei
        vx = (i-1)*dv;
        thx = (i-1)*dv;

        mx = Supersonic('Mnu','v',vx*(180.0/pi),'g',gamma);
        ax = Supersonic('mu','M',mx)*(pi/180.0);

        TriAngle = 0.5*vx+abs(tht);        
        ChordLength = sqrt(2 * ((rrre * rt)^2)*(1 - cos(vx)));

        xr(i) = xr(1) + (ChordLength*cos(TriAngle));
        yr(i) = yr(1) + (ChordLength*sin(TriAngle));

        A = [1 -tan(thxc); 1 -tan(thx-ax)];
        B = [yc(i-1)-tan(thxc)*xc(i-1); yr(i)-tan(thx-ax)*xr(i)];
        solution = A\B;
        yc(i) = solution(1);
        xc(i) = solution(2);

        thxc = thx;
    end

    dv = (ve-vei)/(ne-1.0);
    for i = 1:ne-1
        thxc = thxc - dv;
        vx = vx + dv;

        mx = Supersonic('Mnu','v',vx*(180.0/pi),'g',gamma);
        ax = Supersonic('mu','M',mx)*(pi/180.0);

        A = [1 -tan(thx); 1 -tan(thxc+ax)];
        B = [yr(nei+i-1)-tan(thx)*xr(nei+i-1); yc(nei)-(tan(thxc + ax)*xc(nei))];
        solution = A\B;
        yr(nei+i) = solution(1);
        xr(nei+i) = solution(2);
        thx = thxc;
    end

    rExit = yr(length(yr))
    for i = 1:length(yr)
        yr(i) = rExit - yr(i);
    end
    for i = 1:length(yc)
        yc(i) = rExit - yc(i);
    end

    plot(xr, yr, xc, yc)
end
