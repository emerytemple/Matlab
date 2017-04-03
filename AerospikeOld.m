function[] = Aerospike()
clear all
clc

format long

[leax, leay] = LEA(15,90,1.4,0.1, 50);
[ae2x, ae2y] = AE2();
plot(leax, leay, ae2x, ae2y)
end

function [x,y] = LEA(PR,theta,gamma, throat_len, nr) % linear external (approximate) contour
% Angelino, G. - Approximate Method for Plug Nozzle Design (AIAA)

    x = zeros(nr);
    y = zeros(nr);
    
    mach_exit = Supersonic('Mpr','PR',PR,'g',gamma);
    dm = (mach_exit-1.0)/(nr-1.0);

    % Po = 101325*PR

    for i = 1:nr
        M = 1.0+((i-1)*dm);
        mu = Supersonic('mu','M',M);
        nu = Supersonic('nu','M',M,'g',gamma);

        alpha = (mu+theta-nu)*(pi/180.0);
        xsi = M*Supersonic('A/A*','M',M,'g',gamma);

        % redimensionalize
        L = xsi*throat_len;
        x(i) = L*cos(alpha);
        y(i) = -L*sin(alpha);
    end
end

function [] = LE() % linear external contour
% BDenton axisymmetric 1

Mexit = 3.2
Gamma = 1.2424
rThroat = 1.0
Beta = 1.0
DeltaVAeroD = 0.1

NumChar = 1;
PointNum = 0;
ii = 1;

MaContinue = 0;
while MaContinue == 0
    
    PointNum = PointNum + (NumChar + 1);
    for ii = (PointNum - NumChar):PointNum
        if ii == (PointNum - NumChar)
            Position(ii) = 1;
        elseif ii == PointNum
            Position(ii) = 3;
        else
            Position(ii) = 2;
        end
    end

    for ii = (PointNum - NumChar):PointNum
        if Position(ii) == 1
            Theta(ii) = NumChar * DeltaVAeroD;
            x(ii) = (Beta * rThroat) * sin(Theta(ii));
            r(ii) = rThroat + ((Beta * rThroat) *(1 - cos(Theta(ii))));
            v(ii) = Theta(ii);
            vRad = v(ii);
            Ma(ii) = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
            Alpha(ii) = Supersonic('mu','M',Ma(ii))*(pi/180.0);
        elseif Position(ii) == 3
            Theta(ii) = 0;
            r(ii) = 0;
            Thetatemp = Theta((ii-1));
            rtemp = r((ii-1));
            xtemp = x((ii-1));
            Alphatemp = Alpha((ii-1));
            Matemp = Ma((ii-1));
            vtemp = v((ii-1));
            DeltaTheta = 1.0;
            ThetaLast = 100;
            while DeltaTheta >= 1e-10
                x(ii) = ((rtemp - (tan(Thetatemp -...
                Alphatemp) * xtemp)) /...
                (-tan(Thetatemp - Alphatemp))) -...
                r(ii);
                v(ii) = Thetatemp + vtemp - Theta(ii) + ...
                ((1/(sqrt((Matemp^2) - 1) -...
                (1/tan(Thetatemp)))) * ((r(ii) -...
                rtemp) / rtemp));
                vRad = v(ii);
                Ma(ii) = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
                Alpha(ii) = Supersonic('mu','M',Ma(ii))*(pi/180.0);
                DeltaTheta = abs(ThetaLast - Theta(ii));
                if DeltaTheta > 1e-10
                    ThetaLast = Theta(ii);
                    Thetatemp = (Thetatemp + Theta(ii)) / 2;
                    Alphatemp = (Alphatemp + Alpha(ii)) / 2;
                    Matemp = (Matemp + Ma(ii)) / 2;
                    rtemp = (rtemp + r(ii)) / 2;
                    xtemp = (xtemp + x(ii)) / 2;
                    vtemp = (vtemp + v(ii)) / 2;
                end
            end
            if Ma(ii) >= Mexit
                MaContinue = 1;
            else
                MaContinue = 0;
            end
            NumChar = NumChar + 1;
        else
            ThetatempRight = Theta((ii-1));
            rtempRight = r((ii-1));
            xtempRight = x((ii-1));
            AlphatempRight = Alpha((ii-1));
            MatempRight = Ma((ii-1));
            vtempRight = v((ii-1));
            ThetatempLeft = Theta((ii-NumChar));
            rtempLeft = r((ii-NumChar));
            xtempLeft = x((ii-NumChar));
            AlphatempLeft = Alpha((ii-NumChar));
            MatempLeft = Ma((ii-NumChar));
            vtempLeft = v((ii-NumChar));
            DeltaTheta = 1.0;
            ThetaLast = 100;
            while DeltaTheta >= 1e-10
                a = tan(ThetatempRight - AlphatempRight);
                b = tan(ThetatempLeft + AlphatempLeft);
                c = rtempRight - (a * xtempRight);
                d = rtempLeft - (b * xtempLeft);
                A = [1 -a; 1 -b];
                B = [c; d];
                solution = A\B;
                r(ii) = solution(1);
                x(ii) = solution(2);
                c = (ThetatempRight + vtempRight) + ((1 /...
                ((sqrt(MatempRight^2 - 1)) -...
                (1/tan(ThetatempRight)))) * ((r(ii) -...
                rtempRight) / rtempRight));
                if ThetatempLeft == 0.0
                    d = (2 * ThetatempLeft) - vtempLeft;
                    A = [1 1; 2 -1];
                else
                    d = (ThetatempLeft - vtempLeft) - ((1 /...
                    ((sqrt(MatempLeft^2 - 1)) +...
                    (1/tan(ThetatempLeft)))) * ((r(ii) -...
                    rtempLeft) / rtempLeft));
                    A = [1 1; 1 -1];
                end
                B = [c; d];
                solution = A\B;
                Theta(ii) = solution(1);
                v(ii) = solution(2);
                vRad = v(ii);
                Ma(ii) = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
                Alpha(ii) = asin(1/Ma(ii));
                DeltaTheta = abs(ThetaLast - Theta(ii));
                if DeltaTheta > 1e-10
                    ThetaLast = Theta(ii);
                    ThetatempRight = (ThetatempRight +...
                    Theta(ii)) / 2;
                    AlphatempRight = (AlphatempRight +...
                    Alpha(ii)) / 2;
                    MatempRight = (MatempRight + Ma(ii)) / 2;
                    rtempRight = (rtempRight + r(ii)) / 2;
                    xtempRight = (xtempRight + x(ii)) / 2;
                    vtempRight = (vtempRight + v(ii)) / 2;
                    ThetatempLeft = (ThetatempLeft +...
                    Theta(ii)) / 2;
                    AlphatempLeft = (AlphatempLeft +...
                    Alpha(ii)) / 2;
                    MatempLeft = (MatempLeft + Ma(ii)) / 2;
                    rtempLeft = (rtempLeft + r(ii)) / 2;
                    xtempLeft = (xtempLeft + x(ii)) / 2;
                    vtempLeft = (vtempLeft + v(ii)) / 2;
                end
            end
        end
    end
end

ii = PointNum + 1;
Theta(ii) = Theta((ii-(NumChar-1)));
a = tan(Theta((ii-(NumChar-1))) + Alpha((ii-(NumChar-1))));
b = (Theta((ii-NumChar)) + Theta((ii-(NumChar-1)))) / 2;
c = r((ii-(NumChar-1))) - (a * x((ii-(NumChar-1))));
d = r((ii-NumChar)) - (b * x((ii-NumChar)));
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
r(ii) = solution(1);
x(ii) = solution(2);

jj = 1;
for jj = (ii + 1):(PointNum+(NumChar-1))
    Theta(jj) = Theta((jj-(NumChar-1)));
    a = tan(Theta((jj-(NumChar-1))) + Alpha((jj-(NumChar-1))));
    b = (Theta((jj-1)) + Theta((jj-(NumChar-1)))) / 2;
    c = r((jj-(NumChar-1))) - (a * x((jj-(NumChar-1))));
    d = r((jj-1)) - (b * x((jj-1)));
    A = [1 -a; 1 -b];
    B = [c; d];
    solution = A\B;
    r(jj) = solution(1);
    x(jj) = solution(2);
end

jj = 1;
Char = 1;
for ii = 1:1:(NumChar-1)
    xcontour(ii) = x(jj);
    rcontour(ii) = r(jj);
    jj = jj + (Char + 1);
    Char = Char + 1;
end

jj = PointNum + 1;
for ii = NumChar:((NumChar-1)+(NumChar-1))
    xcontour(ii) = x(jj);
    rcontour(ii) = r(jj);
    jj = jj + 1;
end

xcontour'
rcontour'

% plot(xcontour,rcontour)
end

function [] = LIE() % linear internal-external contour
% bdenton IEaerospike (axisymmetric???)

Mexit = 3.2
Gamma = 1.2424
rThroat = 1.0
Beta = 1.0
DeltaVAeroD = 0.1
Truncate = 0

vmax = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mexit^2) - 1))) - atan((sqrt((Mexit^2) - 1)))
Mei = 2.0
vRegionOne = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mei^2) - 1))) - atan((sqrt((Mei^2) - 1)))
ThetaAeroThroat = vmax - (2 * vRegionOne)
xAeroStreamCowl(1,1) = 0.0
rAeroStreamCowl(1,1) = 0.0
ThetaAeroStreamCowl(1,1) = ThetaAeroThroat
vAeroStreamCowl(1,1) = 0.0
MaAeroStreamCowl(1,1) = 1.0
AlphaAeroStreamCowl(1,1) = asin(1/MaAeroStreamCowl(1,1))
ThetaAeroStream(1,1) = ThetaAeroThroat
vAeroStream(1,1) = 0.0
MaAeroStream(1,1) = 1.0
AlphaAeroStream(1,1) = asin(1/MaAeroStream(1,1))
if ThetaAeroThroat < 0.0
    xAeroStream(1,1) = rThroat * sin(ThetaAeroThroat)
    rAeroStream(1,1) = rThroat * cos(ThetaAeroThroat)
elseif ThetaAeroThroat > 0.0
    xAeroStream(1,1) = -rThroat * sin(ThetaAeroThroat)
    rAeroStream(1,1) = rThroat * cos(ThetaAeroThroat)
else
    xAeroStream(1,1) = 0.0
    rAeroStream(1,1) = rThroat
end
if ThetaAeroThroat < 0.0
    xCenter = (rThroat + (Beta * rThroat)) * sin(ThetaAeroThroat)
    rCenter = (rThroat + (Beta * rThroat)) * cos(ThetaAeroThroat)
elseif ThetaAeroThroat > 0.0
    xCenter = -(rThroat + (Beta * rThroat)) * sin(ThetaAeroThroat)
    rCenter = (rThroat + (Beta * rThroat)) * cos(ThetaAeroThroat)
else
    xCenter = 0.0
    rCenter = rThroat + (Beta * rThroat)
end
ii = 1
jj = 1
vRegionOneCheck = 1
while vRegionOneCheck == 1
    ii = ii + 1
    jj = jj + 1
    ThetaAeroStream(ii,1) = ThetaAeroStream((ii-1),1) + DeltaVAeroD
    vAeroStream(ii,1) = vAeroStream((ii-1),1) + DeltaVAeroD
    if vAeroStream(ii,1) > vRegionOne
        vAeroStream(ii,1) = vRegionOne
        ThetaAeroStream(ii,1) = ThetaAeroThroat + vRegionOne
    end
    vRad = vAeroStream(ii,1)
    MaAeroStream(ii,1) = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
    AlphaAeroStream(ii,1) = asin(1/MaAeroStream(ii,1))
    if ThetaAeroThroat < 0.0
        TriAngle = ((pi/2) - ThetaAeroThroat) - ((pi -...
        vAeroStream(ii,1)) / 2)
    elseif ThetaAeroThroat > 0.0
        TriAngle = ((pi/2) + ThetaAeroThroat) - ((pi -...
        vAeroStream(ii,1)) / 2)
    else
        TriAngle = (pi/2) - ((pi - vAeroStream(ii,1)) / 2)
    end
    ChordLength = sqrt(2 * ((Beta * rThroat)^2) *...
    (1 - cos(vAeroStream(ii,1))))
    DeltaR = ChordLength * sin(TriAngle)
    DeltaX = ChordLength * cos(TriAngle)
    xAeroStream(ii,1) = xAeroStream(1,1) + DeltaX
    rAeroStream(ii,1) = rAeroStream(1,1) + DeltaR
    LineSlope = ThetaAeroStream(ii,1) - AlphaAeroStream(ii,1)
    ThetaLast = ThetaAeroStreamCowl((jj-1),1)
    xLast = xAeroStreamCowl((jj-1),1)
    rLast = rAeroStreamCowl((jj-1),1)
    a = -tan(ThetaLast)
    b = -tan(LineSlope)
    c = rLast - tan(ThetaLast) * xLast
    d = rAeroStream(ii,1) - tan(LineSlope)*xAeroStream(ii,1)
    A = [1 a; 1 b]
    B = [c; d]
    solution = A\B
    rAeroStreamCowl(jj,1) = solution(1,1)
    xAeroStreamCowl(jj,1) = solution(2,1)
    ThetaAeroStreamCowl(jj,1) = ThetaAeroStream(ii,1)
    vAeroStreamCowl(jj,1) = vAeroStream(ii,1)
    MaAeroStreamCowl(jj,1) = MaAeroStream(ii,1)
    AlphaAeroStreamCowl(jj,1) = AlphaAeroStream(ii,1)
    if vAeroStream(ii) >= vRegionOne
        vRegionOneCheck = 0
    else
        vRegionOneCheck = 1
    end
end
MaContinue = 1
DeltaVAeroTemp = 0.0
while MaContinue == 1
    ii = ii + 1;
    ThetaAeroStreamCowl(jj,1) = ThetaAeroStreamCowl(jj,1) - DeltaVAeroTemp
    vAeroStreamCowl(jj,1) = vAeroStreamCowl(jj,1) + DeltaVAeroTemp
    vRad = vAeroStreamCowl(jj,1)
    MaAeroStreamCowl(jj,1) = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
    AlphaAeroStreamCowl(jj,1) = asin(1/MaAeroStreamCowl(jj,1))
    LineSlope = ThetaAeroStreamCowl(jj,1) + AlphaAeroStreamCowl(jj,1)
    rStart = rAeroStream((ii-1),1)
    rIntercept = rAeroStreamCowl(jj,1) -...
    (tan(LineSlope) * xAeroStreamCowl(jj,1))
    ThetaLast = ThetaAeroStream((ii-1),1)
    xLast = xAeroStream((ii-1),1)
    rLast = rAeroStream((ii-1),1)
    a = -tan(ThetaLast)
    b = -tan(LineSlope)
    c = rLast - tan(ThetaLast) * xLast
    d = rIntercept
    A = [1 a; 1 b]
    B = [c; d]
    solution = A\B;
    rAeroStream(ii,1) = solution(1,1)
    xAeroStream(ii,1) = solution(2,1)
    ThetaAeroStream(ii,1) = ThetaAeroStreamCowl(jj,1)
    vAeroStream(ii,1) = vAeroStreamCowl(jj,1)
    MaAeroStream(ii,1) = MaAeroStreamCowl(jj,1)
    AlphaAeroStream(ii,1) = AlphaAeroStreamCowl(jj,1)
    DeltaVAeroTemp = DeltaVAeroD
    if vAeroStreamCowl(jj,1) >= vmax
        MaContinue = 0
    else
        MaContinue = 1
    end
end
if Truncate == 1
    zz = 1;
    PercentLength = 0.0; % Initializes loop
    idealLength = xAeroStream(ii,1);
    while PercentLength <= (Percent/100)
        xAeroStreamTrunc(zz,1) = xAeroStream(zz,1);
        rAeroStreamTrunc(zz,1) = rAeroStream(zz,1);
        ThetaAeroStreamTrunc(zz,1) = ThetaAeroStream(zz,1);
        vAeroStreamTrunc(zz,1) = vAeroStream(zz,1);
        MaAeroStreamTrunc(zz,1) = MaAeroStream(zz,1);
        AlphaAeroStreamTrunc(zz,1) = AlphaAeroStream(zz,1);
        zz = zz + 1;
        PercentLength = abs(xAeroStreamTrunc((zz-1),1))/idealLength;
    end
end
rExit = rAeroStream(ii,1)
for ll = 1:1:(ii)
    rAeroStream(ll,1) = rExit - rAeroStream(ll,1)
end
for ll = 1:1:(jj)
    rAeroStreamCowl(ll,1) = rExit - rAeroStreamCowl(ll,1)
end
rCenter = rExit - rCenter
if Truncate == 1
    for jj=1:1:(zz-1)
        rAeroStreamTrunc(jj,1) = rExit - rAeroStreamTrunc(jj,1); 
    end
end

xAeroStream
rAeroStream
xAeroStreamCowl
rAeroStreamCowl
plot(xAeroStream, rAeroStream, xAeroStreamCowl, rAeroStreamCowl)

end

function [] = AE1() % axi-symmetric external contour 1
% bdenton dentonaerospike

Mexit = 3.2
Gamma = 1.2424
rThroat = 1.0
Beta = 1.0
DeltaVAeroD = 0.1
Truncate = 0

vmax = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mexit^2) - 1))) - atan((sqrt((Mexit^2) - 1)))

xAeroStreamContour(1,1) = -rThroat * cos((pi/2) - vmax)
rAeroStreamContour(1,1) = rThroat * sin((pi/2) - vmax)
ThetaAeroStreamContour(1,1) = vmax
vAeroStreamContour(1,1) = 0.0
MaAeroStreamContour(1,1) = 1.0
AlphaAeroStreamContour(1,1) = asin(1/MaAeroStreamContour(1,1))

xAeroExpansion = 0.0
rAeroExpansion = 0.0
ThetaAeroExpansion = vmax
DeltaVAeroD
vAeroExpansion = 0.0
MaAeroExpansion = 1.0
AlphaAeroExpansion = asin(1/MaAeroExpansion)
MaContinue = 1
NumCharUsed = 1
ii = 2
while MaContinue == 1
    ThetaAeroExpansion = ThetaAeroExpansion - DeltaVAeroD
    if ThetaAeroExpansion < 0.0
        ThetaAeroExpansion = 0.0
    else
        ThetaAeroExpansion = ThetaAeroExpansion
    end
    vAeroExpansion = (NumCharUsed * DeltaVAeroD)
    if vAeroExpansion > vmax
        vAeroExpansion = vmax
    else
        vAeroExpansion = vAeroExpansion
    end
    vRad = vAeroExpansion
    MaAeroExpansion = Supersonic('Mnu','v',vRad*(180.0/pi),'g',Gamma);
    AlphaAeroExpansion = asin(1/MaAeroExpansion)
    LineSlope = AlphaAeroExpansion + ThetaAeroExpansion
    rIntercept = 0.0
    ThetaLast = ThetaAeroStreamContour((ii-1),1)
    xLast = xAeroStreamContour((ii-1),1)
    rLast = rAeroStreamContour((ii-1),1)
    a = -tan(ThetaLast)
    b = -tan(LineSlope)
    c = rLast - tan(ThetaLast) * xLast
    d = rIntercept
    A = [1 a; 1 b]
    B = [c; d]
    solution = A\B
    rAeroStreamContour(ii,1) = solution(1,1)
    xAeroStreamContour(ii,1) = solution(2,1)
    ThetaAeroStreamContour(ii,1) = ThetaAeroExpansion
    vAeroStreamContour(ii,1) = vAeroExpansion
    MaAeroStreamContour(ii,1) = MaAeroExpansion
    AlphaAeroStreamContour(ii,1) = asin(1/MaAeroExpansion)
    ii = ii + 1
    NumCharUsed = NumCharUsed + 1
    if ThetaAeroExpansion <= 0.0
        MaContinue = 0
    else
        MaContinue = 1
    end
end

xAeroStreamContour
rAeroStreamContour

if Truncate == 1
    zz = 1;
    PercentLength = 0.0;
    idealLength = xAeroStreamContour((ii-1),1);
    while PercentLength <= (Percent/100)
        xAeroStreamContourTrunc(zz,1) = xAeroStreamContour(zz,1);
        rAeroStreamContourTrunc(zz,1) = rAeroStreamContour(zz,1);
        ThetaAeroStreamContourTrunc(zz,1) = ThetaAeroStreamContour(zz,1);
        vAeroStreamContourTrunc(zz,1) = vAeroStreamContour(zz,1);
        MaAeroStreamContourTrunc(zz,1) = MaAeroStreamContour(zz,1);
        AlphaAeroStreamContourTrunc(zz,1) = AlphaAeroStreamContour(zz,1);
        zz = zz + 1;
        PercentLength = abs(xAeroStreamContourTrunc((zz-1),1))/idealLength;
    end
end
xShift = abs(xAeroStreamContour(1,1));
rExit = rAeroStreamContour((ii-1),1);
xAeroExpansion = xAeroExpansion + xShift;
rAeroExpansion = rAeroExpansion + rAeroStreamContour((ii-1),1);
for jj=1:1:(ii-1)
    xAeroStreamContour(jj,1) = xAeroStreamContour(jj,1) + xShift;
    rAeroStreamContour(jj,1) = rExit - rAeroStreamContour(jj,1);
end

if Truncate == 1
    for jj=1:1:(zz-1)
        xAeroStreamContourTrunc(jj,1) = xAeroStreamContourTrunc(jj,1)...
        + xShift;
        rAeroStreamContourTrunc(jj,1) = rExit -...
        rAeroStreamContourTrunc(jj,1);
    end
end

    xAeroStreamContour
    rAeroStreamContour
    plot(xAeroStreamContour,rAeroStreamContour)
end

function [x, y] = AE2() % axi-symmetric external contour 2
% lee & thompson

xp = 12.0
n = 50
r = 1718.0
te = 3000.0
papc = 0.01
g = 32.17405
gama = 1.2

rm = Supersonic('Mar','x',xp,'g',gama)
ve = Supersonic('nu','M',rm,'g',gama)*(pi/180.0)

delta = (pi/2.0)-ve
sumcg(1) = ((2.0/(gama+1.0))^(gama/(gama-1.0)))*(gama+1.0)*sin(delta)
ht = (xp-sqrt(xp*(xp-sin(delta))))/(xp*sin(delta))
a1 = (gama+1.0)/(2.0*(gama-1.0))
b1 = sqrt(1.0+0.5*(gama-1.0)*rm*rm)
cfo = gama*rm*((2.0/(gama+1.0))^a1)/b1
pcpt = (0.5*(gama+1.0))^(gama/(gama-1.0))
vt = 1.0+0.5*(gama-1.0)*rm*rm
vt = gama*r*te*vt/(0.5*(gama+1.0))
vt = vt*g
vt = sqrt(vt)
spim = (1.0-pcpt*papc)/gama
spim = 1.0+spim
sumim(1) = vt*sin(delta)*spim/g
vacim = 1.0+1.0/gama
sumva(1) = vt*sin(delta)*vacim/g

xn = n
drm = (rm-1.0)/xn
xm(1)=1.0
rxre(1) = 1.0-ht*sin(delta)
xxre(1) = (-ht)*cos(delta)
a3 = (-gama)/(gama-1.0)
pxpc(1) = (1.0+0.5*(gama-1.0)*xm(1)*xm(1))^a3

i=2
while i-n-1 <= 0
        xm(i) = xm(i-1) + drm
        pro = pxpc(i-1)
        rxo = rxre(i-1)
    a = sqrt((gama-1.0)*(xm(i)*xm(i)-1.0)/(gama+1.0))
    b = sqrt((gama+1.0)/(gama-1.0))
    c = sqrt(xm(i)*xm(i)-1.0)
    c = atan(c)
    vx = b*atan(a)-c
    y = 1.0/xm(i)
    ux = atan(y/sqrt(1.0-y*y))
    a2 = (gama+1.0)/(2.0*(gama-1.0))
    b2 = (2.0/(gama+1.0))*(1.0+0.5*(gama-1.0)*xm(i)*xm(i))
    rxre(i) = 1.0-(b2^a2)*sin(ve-vx+ux)/xp
    rxre(i) = sqrt(rxre(i))
    xxre(i) = (1.0-rxre(i))*cos(ve-vx+ux)/sin(ve-vx+ux)
    a3 = (-gama)/(gama-1.0)
    pxpc(i) = (1.0+0.5*(gama-1.0)*xm(i)*xm(i))^a3
    sumcg(i) = sumcg(i-1)+0.5*xp*(pro+pxpc(i))*(rxo*rxo-rxre(i)*rxre(i))
    c0 = pcpt*vt*xp/(g*gama)
    sumim(i) = sumim(i-1)+0.5*c0*(pro+pxpc(i)-2.0*papc)*(rxo*rxo-rxre(i)*rxre(i))
    sumva(i) = sumva(i-1)+0.5*c0*(pro+pxpc(i))*(rxo*rxo-rxre(i)*rxre(i))
    i = i+1
end

        disp('external')
		disp('i	xm	rx/re	xx/re	px/pc	cfvac	impulse vac. impulse')
		for i = 1:length(xm)
            aaa(i,1) = xm(i);
            aaa(i,2) = real(rxre(i));
            aaa(i,3) = real(xxre(i));
            aaa(i,4) = pxpc(i);
            aaa(i,5) = sumcg(i);
            aaa(i,6) = sumim(i);
            aaa(i,7) = sumva(i);
        end
        aaa

x = real(xxre)
y = real(rxre)
end

function [] = AIE() % axi-symmetric internal-external contour
        % lee & thompson

		n1 = 50;
		n2 = 50;

		peipc = 0.049478861660970451;
		xp = 12.0;
		rrre = 0.3;
		pht = pi/2.0;
		r = 1718.0;
		te = 3000.0;
		papc = 0.0119;
		g = 32.17405;
		gamma = 1.1343;

		gp1 = gamma+1.0;
		gm1 = gamma-1.0;
		exp = gp1/(2.0*gm1);
		c1 = 2.0/gp1;
		c2 = gm1/gp1;
		c3 = gp1/gm1;
		c4 = (1.0-gamma)/gamma;
		c5 = 2.0/gm1;
		c6 = gamma/gm1;

        rm = Supersonic('Mar','x',xp,'g',gamma);
		ve = Supersonic('nu','M',rm,'g',gamma)*(pi/180.0);

		rmei = sqrt(c5*((peipc^c4)-1.0));
		vei = Supersonic('nu','M',rmei,'g',gamma)*(pi/180.0);

		thei = vei - ve;
		phei = thei + atan((1/rmei)/sqrt(1-((1/rmei)*(1/rmei))));

		rpre = sqrt(1-(c1*(1+(0.5*gm1*rmei*rmei)))*sin(phei)/xp);
		xpre = (rpre-1)*cos(phei)/sin(phei);
		drm = (rmei-1.0)/n1;

		% internal portion of the nozzle
% 		xm(1) = 1;
		for k = 1:n1+1
            if k == 1
                xm(1) = 1;
            end
			vx = Supersonic('nu','M',xm(k),'g',gamma)*(pi/180.0);
			bx = pht-(pi/2)-vx+abs(thei);
			x1re = 2.0*rrre*sin(0.5*bx);
			psi = pi-pht+vx-0.5*(pi-bx);
			rx1re(k) = rpre+x1re*sin(psi);
			xx1re(k) = xpre-x1re*cos(psi);

			if k == 1
				rx2re(1) = sqrt(rx1re(1)*rx1re(1)+sin(pht)/xp);
				xx2re(1) = xx1re(1)+(rx2re(1)-rx1re(1))*cos(pht)/sin(pht); % use atan
            else
				phx = 2.0*vei-ve-vx+atan(1.0/(xm(k)*sqrt(1.0-(1.0/xm(k))^2.0)));
				a2 = c1*(1.0+0.5*gm1*xm(k)*xm(k));
				rx2re(k) = sqrt(rx1re(k)*rx1re(k)+(a2^(0.5*c3))*sin(phx)/xp);
				xx2re(k) = xx1re(k)+(rx2re(k)-rx1re(k))*cos(phx)/sin(phx);
            end
			pxpc(k) = (1.0+0.5*gm1*xm(k)*xm(k))^(-gamma/gm1);
			if k < n1+1
				xm(k+1) = xm(k) + drm;
            end
        end

        disp('internal')
		disp('i	xm	rx1re	xx1re	rx2re	xx2re	pxpc')
		for i = 1:n1+1
            aaa(i,1) = xm(i);
            aaa(i,2) = rx1re(i);
            aaa(i,3) = xx1re(i);
            aaa(i,4) = rx2re(i);
            aaa(i,5) = xx2re(i);
            aaa(i,6) = pxpc(i);
        end
        aaa
        
		% external portion of the nozzle
		xm2(1) = xm(k);
		pxpc2(1) = pxpc(k);

		ux = atan(1.0/(xm2(1)*sqrt(1.0-(1.0/xm2(1))^2.0)));
		vx = Supersonic('nu','M',xm2(1),'g',gamma)*(pi/180.0);
		rxre(1) = sqrt(1-((c1*(1+0.5*gm1*xm2(1)*xm2(1)))^exp)*sin(ve-vx+ux)/xp);
		xxre(1) = (1-rxre(1))*cos(ve-vx+ux)/sin(ve-vx+ux);

		ch2 = sqrt((0.5*gp1*xm2(1)*xm2(1))/(1.0+0.5*gm1*xm2(1)*xm2(1)));
		sumcg(1) = gamma*(c1^c6)*ch2*cos(thei)+xp*pxpc2(1)*(1-rxre(1)*rxre(1));

		vt = sqrt(gamma*g*r*(te*(1+0.5*gm1*rm*rm))/(0.5*gp1));
		vc = 1+0.5*gm1*xm2(1)*xm2(1);
		vq = sqrt(gamma*r*g*(te*(1+0.5*gm1*rm*rm))/vc);

		a = (1+0.5*gm1*xm2(1)*xm2(1))^(-c6);
		sumim(1) = vq*(cos(thei)+(1-papc*a*sin(phei)/(xm2(1)*xm2(1)))/gamma)/g;
		sumva(1) = vq*(cos(thei)+1/gamma)/g;

		drm = (rm-xm2(1))/n2;
		xm2(2) = xm2(1)+drm;
		pro = pxpc2(1);
		rxo = rxre(1);

		for k =1:n2
			ux = atan(1/(xm2(k+1)*sqrt(1-(1/xm2(k+1))^2)));
			vx = Supersonic('nu','M',xm2(k+1),'g',gamma)*(pi/180.0);

			rxre(k+1) = c1*(1+0.5*gm1*xm2(k+1)*xm2(k+1));
			rxre(k+1) = rxre(k+1)^(gp1*0.5/gm1);
			rxre(k+1) = sqrt(abs(1-rxre(k+1)*sin(ve-vx+ux)/xp)); % abs is temporary

			xxre(k+1) = (1-rxre(k+1))*cos(ve-vx+ux)/sin(ve-vx+ux);

			pxpc2(k+1) = (1+0.5*gm1*xm2(k+1)*xm2(k+1))^(-c6);

			sumcg(k+1) = sumcg(k)+0.5*xp*(pro+pxpc2(k+1))*(rxo*rxo-rxre(k+1)*rxre(k+1));

			a = c6;
			a = (0.5*gp1)^a;
			a = a*vt/(gamma*g);
			b = 0.5*a*xp;
			sumim(k+1) = sumim(k)+b*(pro+pxpc2(k+1)-2*papc)*(rxo*rxo-rxre(k+1)*rxre(k+1));
			sumva(k+1) = sumva(k)+b*(pro+pxpc2(k+1))*(rxo*rxo-rxre(k+1)*rxre(k+1));

			xm2(k+2) = xm2(k+1) + drm;

			pro = pxpc2(k+1);
			rxo = rxre(k+1);
        end
        
		disp('external')
		disp('i	xm2	rxre	xxre	pxpc2	sumcg   sumim  sumva')
		for i = 1:n1+1
            aaa(i,1) = xm2(i);
            aaa(i,2) = rxre(i);
            aaa(i,3) = xxre(i);
            aaa(i,4) = pxpc2(i);
            aaa(i,5) = sumcg(i);
            aaa(i,6) = sumim(i);
            aaa(i,7) = sumva(i);
        end
        aaa
end




















