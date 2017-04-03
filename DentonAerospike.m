%**************************************************************************
%* AxiSymmetric Aerospike Nozzle Design MOC
%* Brandon Denton
%* RIT Graduate Student
%**************************************************************************
%**************************************************************************
%* Dictionary of Variables
%**************************************************************************
%* 1. a, b, c, d -------------------- Variables that temporarily hold
%* numbers used to evaluate the point
%* that satisifies the Stream Function
%* 2. A, B -------------------------- Matrices that hold the equations to
%* solve for the position of the point
%* that satisfies the Stream Function
%* 3. AlphaAeroExpansion ------------ The Mach angle of the fluid at the
%* expansion point
%* 4. AlphaAeroStreamContour() ------ Array that holds the Mach angle of
%* the points on the contour of the
%* Aerospike nozzle
%* 5. AlphaAeroStreamContourTrunc() - Array that holds the Mach angle of
%* the points on the contour of the
%* truncated Aerospike nozzle
%* 6. DeltaVAeroD ------------------- The desired change in Prandtl-Meyer
%* Expansion Angle
%* 7. Gamma ------------------------- Ratio of Specific Heats of the
%* working fluid
%* 8. idealLength ------------------- Variable that holds the length of an
%* idealized (100% Length) Aerospike
%* nozzle
%* 9. ii, jj ------------------------ Loop and position indices
%* 10. LineSlope --------------------- The slope of the line in radians
%* that approximates the Mach Line in
%* the flowfield
%* 11. MaAeroExpansion --------------- Mach Number of the fluid at the
%* expansion point
%* 12. MaAeroStreamContour() --------- Array that holds the Mach Number of
%* the points on the contour of the
%* Aerospike nozzle
%* 13. MaAeroStreamContourTrunc() ---- Array that holds the Mach Number of
%* the points on the contour of the
%* Aerospike nozzle
%* 14. MaContinue -------------------- Control variable for the loop that
%* calculates the Aerospike nozzle's
%* contour
%* 15. Mach -------------------------- Solution variable of the Mach Number
%* returned from subroutine PMtoMA
%* 16. MachG ------------------------- Initial guess value of the Mach
%* Number used in the Subroutine PMtoMA
%* to find the point's Mach Number
%* 17. Mexit ------------------------- Desired Exit Mach Number
%* 18. NumCharUsed ------------------- The number of C+ characteristics
%* used in the calculation of the
%* Aerospike's contour
%* 19. PercentLength ----------------- Control variable that is calculated
%* to make sure that the Aerospike
%* nozzle is truncate appropriately
%* 20. PMtoMA ------------------------ Subroutine used to find the Mach
%* Number of the point in the flowfield
%* 21. rAeroExpansion ---------------- The r-component of the position of
%* the expansion point of the Aerospike
%* nozzle
%* 22. rAeroStreamContour() ---------- Array that holds the r-component of
%* the points on the contour of the
%* Aerospike nozzle
%* 23. rAeroStreamContourTrunc() ----- Array that holds the r-component of
%* the points on the contour of the
%* truncated Aerospike nozzle
%* 24. rExit ------------------------- Variable used to "flip" the
%* Aerospike's contour so that the axis
%* of symmetry of the nozzle is the
%* x-axis
%* 25. rIntercept -------------------- The r-intercept of the line that
%* approximates the Mach Line in the
%* flowfield
%* 26. rLast ------------------------- The variable that holds the
%* r-component of the point that last
%* satisfied the Stream Function
%* 27. rThroat ----------------------- User-defined radius of the throat
%* 28. solution() -------------------- Matrix that holds the solution to
%* the coordinates of the point that
%* satisfies the Stream Function
%* 29. ThetaAeroExpansion ------------ The flow direction of the fluid at
%* the expansion point
%* 30. ThetaAeroStreamContour() ------ Array that holds the flow direction
%* of the point on the contour of the
%* Aerospike nozzle
%* 31. ThetaAeroStreamContourTrunc() - Array that holds the flow direction
%* of the point on the contour of the
%* truncated Aerospike nozzle
%* 32. ThetaLast --------------------- The variable the holds the flow
%* direction of the point that last
%* satisfied the Stream Function
%* 33. Truncate ---------------------- Control variable that controls
%* whether or not the user wanted to
%* truncate the Aerospike nozzle
%* 34. vAeroStreamContour() ---------- Array that holds the Prandtl-Meyer
%* expansion angle of the point on the
%* contour of the Aerospike nozzle
%* 35. vAeroStreamContourTrunc() ----- Array that holds the Prandtl-Meyer
%* expansion angle of the point on the
%* contour of the Aerospike nozzle
%* 36. vmax -------------------------- Prandtl-Meyer Expansion Angle for
%* the desired Mach Number
%* 37. vRad -------------------------- Variable used in subroutine PMtoMA
%* to calculate the Mach number at the
%* point in question of the flow
%* 38. xAeroExpansion ---------------- The x-component of the position of
%* the expansion point of the Aerospike
%* nozzle
%* 39. xAeroStreamContour() ---------- Array that holds the x-component of
%* the points on the contour of the
%* Aerospike nozzle
%* 40. xAeroStreamContourTrunc() ----- Array that holds the x-component of
%* the points on the contour of the
%* truncated Aerospike nozzle
%* 41. xLast ------------------------- The variable that holds the
%* x-component of the point that last
%* satisfied the Stream Function
%* 42. xShift ------------------------ Variable used to shift the Aerospike
%* contour so that the first point is
%* on the r-axis
%* 43. zz ---------------------------- Indice variable
%**************************************************************************
%**************************************************************************
%* Start Program
%**************************************************************************
%**************************************************************************
format long
disp(' ');
% Calculate the Maximum Prandtl-Meyer Expansion Angle
vmax = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mexit^2) - 1))) - atan((sqrt((Mexit^2) - 1)));
% Calculate the Position of the end points of the sonic line
xAeroStreamContour(1,1) = -rThroat * cos((pi/2) - vmax);
rAeroStreamContour(1,1) = rThroat * sin((pi/2) - vmax);
ThetaAeroStreamContour(1,1) = vmax;
vAeroStreamContour(1,1) = 0.0;
MaAeroStreamContour(1,1) = 1.0;
AlphaAeroStreamContour(1,1) = asin(1/MaAeroStreamContour(1,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressExtAero(1,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStreamContour(1,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempExtAero(1,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStreamContour(1,1)) ^ 2))) ^ (-1);
DeltaSExtAero(1,1) = cp * log(TempExtAero(1,1)/Tstag) - Rgas *...
log(PressExtAero(1,1)/Pstag);
end
% Position of the Expansion point
xAeroExpansion = 0.0;
rAeroExpansion = 0.0;
ThetaAeroExpansion = vmax;
vAeroExpansion = 0.0;
MaAeroExpansion = 1.0;
AlphaAeroExpansion = asin(1/MaAeroExpansion);
% Calculate the Aerospike contour by assuming that the C+ Characteristice
% eminating from the expansion point is a straight line and terminates
% when it satisfies the streamline condition (becomes the point on the
% contour)
MaContinue = 1;
NumCharUsed = 1;
ii = 2;
while MaContinue == 1
% Set expansion point conditions
ThetaAeroExpansion = ThetaAeroExpansion - DeltaVAeroD;
if ThetaAeroExpansion < 0.0
ThetaAeroExpansion = 0.0;
else
ThetaAeroExpansion = ThetaAeroExpansion;
end
vAeroExpansion = (NumCharUsed * DeltaVAeroD);
if vAeroExpansion > vmax
vAeroExpansion = vmax;
else
vAeroExpansion = vAeroExpansion;
end
%Calculate the Mach Number at the current point
vRad = vAeroExpansion;
MachG = MaAeroExpansion;
PMtoMA %calls subprogram to find the Mach Number
MaAeroExpansion = Mach;
AlphaAeroExpansion = asin(1/MaAeroExpansion);
LineSlope = AlphaAeroExpansion + ThetaAeroExpansion;
rIntercept = 0.0; % Center of the expansion wave is located at (0,0)
% Calculate the point that satisfies the Stream Function condition at
% the last expansion point on the contour
% Initiate values of the last point on the streamline
ThetaLast = ThetaAeroStreamContour((ii-1),1);
xLast = xAeroStreamContour((ii-1),1);
rLast = rAeroStreamContour((ii-1),1);
% Calculate the position of the point the satisfies the streamline
% condition
a = -tan(ThetaLast);
b = -tan(LineSlope);
c = rLast - tan(ThetaLast) * xLast;
d = rIntercept;
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rAeroStreamContour(ii,1) = solution(1,1);
xAeroStreamContour(ii,1) = solution(2,1);
ThetaAeroStreamContour(ii,1) = ThetaAeroExpansion;
vAeroStreamContour(ii,1) = vAeroExpansion;
MaAeroStreamContour(ii,1) = MaAeroExpansion;
AlphaAeroStreamContour(ii,1) = asin(1/MaAeroExpansion);
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressExtAero(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStreamContour(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempExtAero(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStreamContour(ii,1)) ^ 2))) ^ (-1);
DeltaSExtAero(ii,1) = cp * log(TempExtAero(ii,1)/Tstag) - Rgas *...
log(PressExtAero(ii,1)/Pstag);
end
ii = ii + 1; %Increase the index
NumCharUsed = NumCharUsed + 1;
if ThetaAeroExpansion <= 0.0
MaContinue = 0;
else
MaContinue = 1;
end
end
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressExtAero(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroExpansion) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempExtAero(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroExpansion) ^ 2))) ^ (-1);
DeltaSExtAero(ii,1) = cp * log(TempExtAero(ii,1)/Tstag) - Rgas *...
log(PressExtAero(ii,1)/Pstag);
end
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
% The program will now "flip" the contour so that the axisymmetric line
% aligns with r=0 and the expansion point will occur at (0,rexit). it will
% also shift the nozzle to that the beginning of the contour will exit on
% the x=0 axis
xShift = abs(xAeroStreamContour(1,1));
rExit = rAeroStreamContour((ii-1),1);
xAeroExpansion = xAeroExpansion + xShift;
rAeroExpansion = rAeroExpansion + rAeroStreamContour((ii-1),1);
for jj=1:1:(ii-1)
xAeroStreamContour(jj,1) = xAeroStreamContour(jj,1) + xShift;
rAeroStreamContour(jj,1) = rExit - rAeroStreamContour(jj,1);
end
% This section only calculates if the Aerospike is truncated
if Truncate == 1
for jj=1:1:(zz-1)
xAeroStreamContourTrunc(jj,1) = xAeroStreamContourTrunc(jj,1)...
+ xShift;
rAeroStreamContourTrunc(jj,1) = rExit -...
rAeroStreamContourTrunc(jj,1);
end
end