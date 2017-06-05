function [xc, yc, xr, yr] = linear_internal() %(tht, xpe, rrre, gamma, dv)
    gamma = 1.4;
    xpe = 5.0;
    tht = 0.0;
    rt = 1.0;
    Beta = 0.5;
    dv = 0.001;
    nei = 20;

    me = Supersonic('Mar','x',xpe,'g',gamma);
    ve = Supersonic('nu','M',me,'g',gamma);

    vei = (ve-tht)/2.0;
    mei = Supersonic('Mnu','v',vei,'g',gamma);

    ve = ve*(pi/180.0);
    vei = vei*(pi/180.0);
    tht = tht*(pi/180.0);

    xc(1) = 0.0;
    yc(1) = 0.0;
    
    ThetaAeroStreamCowl(1) = tht;
    vAeroStreamCowl(1) = 0.0;
    MaAeroStreamCowl(1) = 1.0;
    AlphaAeroStreamCowl(1) = Supersonic('mu','M',MaAeroStreamCowl(1))*(pi/180.0);

    ThetaAeroStream(1) = tht;
    vAeroStream(1) = 0.0;
    MaAeroStream(1) = 1.0;
    AlphaAeroStream(1) = Supersonic('mu','M',MaAeroStreamCowl(1))*(pi/180.0);

    if tht < 0.0
        xr(1) = rt * sin(tht);
        yr(1) = rt * cos(tht);

        xCenter = (rt + (Beta * rt)) * sin(tht);
        ycenter = (rt + (Beta * rt)) * cos(tht);
    elseif tht > 0.0
        xr(1) = -rt * sin(tht);
        yr(1) = rt * cos(tht);

        xCenter = -(rt + (Beta * rt)) * sin(tht);
        ycenter = (rt + (Beta * rt)) * cos(tht);
    else
        xr(1) = 0.0;
        yr(1) = rt;

        xCenter = 0.0;
        ycenter = rt + (Beta * rt);
    end

    ii = 1;
    jj = 1;
    veiCheck = 1;

%     while veiCheck == 1
    dv = (mei-1.0)/(nei-1.0);
    for i = 1:nei
        ii = ii + 1;
        jj = jj + 1;

        ThetaAeroStream(ii) = ThetaAeroStream((ii-1)) + dv;
        vAeroStream(ii) = vAeroStream((ii-1)) + dv;
%         if vAeroStream(ii) > vei
%             vAeroStream(ii) = vei;
%             ThetaAeroStream(ii) = tht + vei;
%         end

        MaAeroStream(ii) = Supersonic('Mnu','v',vAeroStream(ii)*(180.0/pi),'g',gamma);
        AlphaAeroStream(ii) = Supersonic('mu','M',MaAeroStream(ii))*(pi/180.0);

        if tht < 0.0
            TriAngle = ((pi/2) - tht) - ((pi -vAeroStream(ii)) / 2);
        elseif tht > 0.0
            TriAngle = ((pi/2) + tht) - ((pi -vAeroStream(ii)) / 2);    
        else
            TriAngle = (pi/2) - ((pi - vAeroStream(ii)) / 2);
        end
        ChordLength = sqrt(2 * ((Beta * rt)^2)*(1 - cos(vAeroStream(ii))));

        xr(ii) = xr(1) + (ChordLength*cos(TriAngle));
        yr(ii) = yr(1) + (ChordLength*sin(TriAngle));

        LineSlope = ThetaAeroStream(ii) - AlphaAeroStream(ii);

        a = -tan(ThetaAeroStreamCowl((jj-1)));
        b = -tan(LineSlope);
        c = yc((jj-1)) - tan(ThetaAeroStreamCowl((jj-1))) * xc((jj-1));
        d = yr(ii) - tan(LineSlope)*xr(ii);
        A = [1 a; 1 b];
        B = [c; d];
        solution = A\B;
        yc(jj) = solution(1);
        xc(jj) = solution(2);
        ThetaAeroStreamCowl(jj) = ThetaAeroStream(ii);
        vAeroStreamCowl(jj) = vAeroStream(ii);
        MaAeroStreamCowl(jj) = MaAeroStream(ii);
        AlphaAeroStreamCowl(jj) = AlphaAeroStream(ii);

%         if vAeroStream(ii) >= vei
%             veiCheck = 0;
%         else
%             veiCheck = 1;
%         end
    end

    MaContinue = 1;
    DeltaVAeroTemp = 0.0;

    while MaContinue == 1
        ii = ii + 1;

        ThetaAeroStreamCowl(jj) = ThetaAeroStreamCowl(jj) - DeltaVAeroTemp;
        vAeroStreamCowl(jj) = vAeroStreamCowl(jj) + DeltaVAeroTemp;

        MaAeroStreamCowl(jj) = Supersonic('Mnu','v',vAeroStreamCowl(jj)*(180.0/pi),'g',gamma);    
        AlphaAeroStreamCowl(jj) = Supersonic('mu','M',MaAeroStreamCowl(jj))*(pi/180.0);

        LineSlope = ThetaAeroStreamCowl(jj) + AlphaAeroStreamCowl(jj);

        a = -tan(ThetaAeroStream((ii-1)));
        b = -tan(LineSlope);
        c = yr((ii-1)) - tan(ThetaAeroStream((ii-1))) * xr((ii-1));
        d = yc(jj) -(tan(LineSlope) * xc(jj));
        A = [1 a; 1 b];
        B = [c; d];
        solution = A\B;
        yr(ii) = solution(1);
        xr(ii) = solution(2);
        ThetaAeroStream(ii) = ThetaAeroStreamCowl(jj);
        vAeroStream(ii) = vAeroStreamCowl(jj);
        MaAeroStream(ii) = MaAeroStreamCowl(jj);
        AlphaAeroStream(ii) = AlphaAeroStreamCowl(jj);

        DeltaVAeroTemp = dv;
        if vAeroStreamCowl(jj) >= ve
            MaContinue = 0;
        else
            MaContinue = 1;
        end
    end

    rExit = yr(ii);
    for ll = 1:1:(ii)
        yr(ll) = rExit - yr(ll);
    end
    for ll = 1:1:(jj)
        yc(ll) = rExit - yc(ll);
    end

    ycenter = rExit - ycenter;

    plot(xr, yr, xc, yc)
end
