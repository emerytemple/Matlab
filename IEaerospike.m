%**************************************************************************
%* Internal-External Expansion Aerospike
%* Brandon Denton
%* RIT Graduate Student
%* January 31, 2007
%**************************************************************************
%**************************************************************************
%* Dictionary of Terms
%**************************************************************************
%* 1. a, b, c, d ------------ Variables that hold temporary values in the
%* calculation of the point that satisfy the
%* Stream Function condition using matrix
%* algebra
%* 2. A, B ------------------ Matrices that hold the temporary values in
%* the calculation of the points that satisfy
%* the Stream Function condition using matrix
%* algebra
%* 3. AlphaAeroStream()------ Array that holds the Mach angle of the
%* nozzle's wall contour
%* 4. AlphaAeroStreamCowl()-- Array that holds the Mach angle of the
%* nozzle's cowl wall contour
%* 5. AlphaAeroStreamTrunc()- Array that holds the Mach angle of the
%* nozzle's truncated wall contour
%* 6. Beta ------------------ Factor of the throat that defines that
%* radius of the arc that defines the first
%* expansion region of the nozzle
%* 7. ChordLength ----------- Length of the chord of the arc between points
%* used to calculate the position of points on
%* the expansion arc of region one
%* 8. DeltaR ---------------- The change in r-component relative to the
%* first expansion arc point of the throat
%* defining the expansion of region one
%* 9. DeltaVAeroD ----------- User-defined change in Prandtl-Meyer
%* expansion angle
%* 10. DeltaVAeroTemp -------- Variable that temporary changes and stores
%* DeltaVAeroD in order to calculate the first
%* C+ characteristic and flow condition in
%* region two
%* 11. DeltaX ---------------- The change in the x-component relative to
%* the first expansion arc point of the throat
%* defining the expansion of region one
%* 12. EndCheck -------------- Variable that controls the loop to make sure
%* the the end guess of the streamline
%* calculation is larger than the actual
%* streamline solution
%* 13. Gamma ----------------- Ratio of Specific Heats of the working fluid
%* 14. idealLength ----------- The total length of the 100% aerospike
%* nozzle
%* 15. ii, jj, zz ------------ Loop and position indices
%* 16. LineSlope ------------- Variable that holds the radian angle of the
%* Mach line accelerating the flow
%* 17. MaAeroStream() -------- Array that holds the Mach number of the
%* nozzle's wall contour
%* 18. MaAeroStreamCowl()----- Array that holds the Mach number of the
%* nozzle's cowl wall contour
%* 19. MaAeroStreamTrunc()---- Array that holds the Mach number of the
%* nozzle's truncated wall contour
%* 20. MaContinue ------------ Variable that controls the loop that
%* calculates the streamline in the second
%* expansion region
%* 21. Mach ------------------ Mach Number solution for the point returned
%* by subroutine PMtoMA
%* 22. MachG ----------------- Initial guess Mach number used by the
%* subroutine PMtoMA to find the points
%* corresponding Mach Number
%* 23. Mei ------------------- Exit Mach number of the internal expansion
%* of the aerospike nozzle
%* 24. Mexit ----------------- Desired exit Mach number of nozzle
%* 25. NumChar --------------- Number of characteristics used in the first
%* expansion region
%* 26. Percent --------------- Variable that holds the User-Defined desired
%* retained length of the aerospike nozzle
%* 27. PercentLength --------- Variable that temporarily holds the
%* calculated length percentage of the
%* truncated nozzle with respect to the ideal
%* (100%) nozzle
%* 28. rAeroStream() --------- Array that holds the r-component of the
%* nozzle's wall contour
%* 29. rAeroStreamCowl()------ Array that holds the r-component of the
%* nozzle's cowl wall contour
%* 30. rAeroStreamTrunc()----- Array that holds the r-component of the
%* nozzle's truncated wall contour
%* 31. rCalc ----------------- Calculated r-component guess value that
%* satisfies the streamline condition
%* 32. rCenter --------------- r-coordinate of the center of the circle
%* defining the entrance region of the internal
%* expansion region of the aerospike nozzle
%* 33. rExit ----------------- Exit radius of the nozzle. Also used in
%* flipping the nozzle so that the axisymmetric
%* line is on the x-axis
%* 34. rLast ----------------- Variable that temporarily holds the
%* r-component of the last calculated point
%* that satisfies the streamline condition
%* 35. rThroat --------------- User-defined radius of the throat
%* 36. solution()------------- Matrix that holds the solution coordinates
%* of the points that satisfy the Stream
%* Function condition
%* 37. ThetaAeroStream()------ Array that holds the values for the flow
%* direction on the nozzle's wall contour
%* 38. ThetaAeroStreamCowl()-- Array that holds the values for the flow
%* direction for the nozzle's cowl contour
%* 39. ThetaAeroStreamTrunc()- Array that holds the flow direction values
%* of the nozzle's truncated wall contour
%* 40. ThetaAeroThroat ------- Angle the sonic line makes with the x-axis
%* 41. ThetaLast ------------- Variable that temporarily holds the flow
%* direction of the last point that satisfied
%* the streamline condition
%* 42. TriAngle -------------- The angle the chord of the expansion arc in
%* region one makes with the x-axis
%* 43. Truncate -------------- Control variable that tells the program
%* whether or not to calculate a truncated
%* version of the aerospike nozzle
%* 44. vmax ------------------ Prandtl-Meyer Expansion angle for the
%* desired exit Mach Number
%* 45. vRad ------------------ Variable that hold the Prandtl-Meyer
%* expansion angle of the point used by the
%* subroutine PMtoMA to find its corresponding
%* Mach number. In Radians
%* 46. vRegionOne ------------ Maximum Pradtl-Meyer Expansion angle allowed
%* in the internal expansion region of the
%* aerospike nozzle
%* 47. vRegionOneCheck ------- Dummy variable that defines the beginning
%* and end of the internal expansion region
%* calculation loop of the aerospike nozzle
%* 48. vAeroStream() --------- Array that holds the values of the expansion
%* angle for the nozzle's wall contour
%* 49. vAeroStreamCowl()----- Array that holds the values of the expansion
%* angle for the nozzle's cowl contour
%* 50. vAeroStreamTrunc()----- Array that holds the flow direction values
%* for the nozzle's truncated wall contour
%* 51. xAeroStream() --------- Array that holds the x-components of the
%* nozzle's wall contour
%* 52. xAeroStreamCowl()------ Array that holds the x-componen fo the
%* nozzle's cowl wall contour
%* 53. xAeroStreamTrunc()----- Array that holds the x-component of the
%* nozzle's truncated wall contour
%* 54. xCenter --------------- x-coordinate of the center of the circle
%* defining the entrance region of the internal
%* expansion region of the aerospike nozzle
%* 55. xLast ----------------- Variable that temporarily holds the
%* x-component of the last point that satisfies
%* the streamline condition
%**************************************************************************
%**************************************************************************
%* START PROGRAM
%**************************************************************************
%**************************************************************************
format long
% Calculate the Maximum Prandtl-Meyer expansion angle for the desired exit
% Mach Number
vmax = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mexit^2) - 1))) - atan((sqrt((Mexit^2) - 1)));
% Calculate the expansion angle for the first expansion region. I have
% chosen 25% of the expansion to take place in this region. For Maximum
% Prandtl-Meyer expansion angles greater than 45 degrees (pi/2 radians) the
% remainder of the expansion angle less pi/2 radians will be used in the
% first expansion region.
Mei = 2.0;
vRegionOne = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((Mei^2) - 1))) - atan((sqrt((Mei^2) - 1)));
% Calculate the flow direction at the throat
ThetaAeroThroat = vmax - (2 * vRegionOne);
% Initialize the condition at both ends of the throat.
% Farther from axisymmetric line end
xAeroStreamCowl(1,1) = 0.0;
rAeroStreamCowl(1,1) = 0.0;
ThetaAeroStreamCowl(1,1) = ThetaAeroThroat;
vAeroStreamCowl(1,1) = 0.0;
MaAeroStreamCowl(1,1) = 1.0;
AlphaAeroStreamCowl(1,1) = asin(1/MaAeroStreamCowl(1,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressAeroStreamCowl(1,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStreamCowl(1,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempAeroStreamCowl(1,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStreamCowl(1,1)) ^ 2))) ^ (-1);
DeltaSAeroStreamCowl(1,1) = cp * log(TempAeroStreamCowl(1,1)/Tstag)...
- Rgas * log(PressAeroStreamCowl(1,1)/Pstag);
end
% Axisymmetric line closer end
ThetaAeroStream(1,1) = ThetaAeroThroat;
vAeroStream(1,1) = 0.0;
MaAeroStream(1,1) = 1.0;
AlphaAeroStream(1,1) = asin(1/MaAeroStream(1,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressAeroStream(1,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStream(1,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempAeroStream(1,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStream(1,1)) ^ 2))) ^ (-1);
DeltaSAeroStream(1,1) = cp * log(TempAeroStream(1,1)/Tstag) - Rgas *...
log(PressAeroStream(1,1)/Pstag);
end
if ThetaAeroThroat < 0.0
xAeroStream(1,1) = rThroat * sin(ThetaAeroThroat);
rAeroStream(1,1) = rThroat * cos(ThetaAeroThroat);
elseif ThetaAeroThroat > 0.0
xAeroStream(1,1) = -rThroat * sin(ThetaAeroThroat);
rAeroStream(1,1) = rThroat * cos(ThetaAeroThroat);
else
xAeroStream(1,1) = 0.0;
rAeroStream(1,1) = rThroat;
end
% Calculate the center of the circle that defines the arc of the expansion
% region
if ThetaAeroThroat < 0.0
xCenter = (rThroat + (Beta * rThroat)) * sin(ThetaAeroThroat);
rCenter = (rThroat + (Beta * rThroat)) * cos(ThetaAeroThroat);
elseif ThetaAeroThroat > 0.0
xCenter = -(rThroat + (Beta * rThroat)) * sin(ThetaAeroThroat);
rCenter = (rThroat + (Beta * rThroat)) * cos(ThetaAeroThroat);
else
xCenter = 0.0;
rCenter = rThroat + (Beta * rThroat);
end
% Initialize loop and indices conditions
ii = 1;
jj = 1;
vRegionOneCheck = 1;
% Loop that calculates point values and streamline points for the cowl wall
% contour
while vRegionOneCheck == 1
ii = ii + 1; %Increase indices every loop iteration
jj = jj + 1;
% Calculate the position of the Region One expansion wall contour
ThetaAeroStream(ii,1) = ThetaAeroStream((ii-1),1) + DeltaVAeroD;
vAeroStream(ii,1) = vAeroStream((ii-1),1) + DeltaVAeroD;
if vAeroStream(ii,1) > vRegionOne
vAeroStream(ii,1) = vRegionOne;
ThetaAeroStream(ii,1) = ThetaAeroThroat + vRegionOne;
end
%Calculate the Mach Number of the current point
vRad = vAeroStream(ii,1);
MachG = 1.0; % Initial Guess of the point's Mach Number
PMtoMA % Calls subroutine to calculate the point's Mach Number
MaAeroStream(ii,1) = Mach;
AlphaAeroStream(ii,1) = asin(1/MaAeroStream(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressAeroStream(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStream(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempAeroStream(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStream(ii,1)) ^ 2))) ^ (-1);
DeltaSAeroStream(ii,1) = cp * log(TempAeroStream(ii,1)/Tstag)...
- Rgas * log(PressAeroStream(ii,1)/Pstag);
end
% Calculate the position of the points defining the entrance region
if ThetaAeroThroat < 0.0
TriAngle = ((pi/2) - ThetaAeroThroat) - ((pi -...
vAeroStream(ii,1)) / 2);
elseif ThetaAeroThroat > 0.0
TriAngle = ((pi/2) + ThetaAeroThroat) - ((pi -...
vAeroStream(ii,1)) / 2);
else
TriAngle = (pi/2) - ((pi - vAeroStream(ii,1)) / 2);
end
ChordLength = sqrt(2 * ((Beta * rThroat)^2) *...
(1 - cos(vAeroStream(ii,1))));
DeltaR = ChordLength * sin(TriAngle);
DeltaX = ChordLength * cos(TriAngle);
xAeroStream(ii,1) = xAeroStream(1,1) + DeltaX;
rAeroStream(ii,1) = rAeroStream(1,1) + DeltaR;
% For the C- Characteristics in Region One
LineSlope = ThetaAeroStream(ii,1) - AlphaAeroStream(ii,1);
% Initiate the values of the last point on the cowl streamline
ThetaLast = ThetaAeroStreamCowl((jj-1),1);
xLast = xAeroStreamCowl((jj-1),1);
rLast = rAeroStreamCowl((jj-1),1);
% Calculate the point that satisfies the Stream Function
a = -tan(ThetaLast);
b = -tan(LineSlope);
c = rLast - tan(ThetaLast) * xLast;
d = rAeroStream(ii,1) - tan(LineSlope)*xAeroStream(ii,1);
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rAeroStreamCowl(jj,1) = solution(1,1);
xAeroStreamCowl(jj,1) = solution(2,1);
ThetaAeroStreamCowl(jj,1) = ThetaAeroStream(ii,1);
vAeroStreamCowl(jj,1) = vAeroStream(ii,1);
MaAeroStreamCowl(jj,1) = MaAeroStream(ii,1);
AlphaAeroStreamCowl(jj,1) = AlphaAeroStream(ii,1);
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressAeroStreamCowl(jj,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStreamCowl(jj,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempAeroStreamCowl(jj,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStreamCowl(jj,1)) ^ 2))) ^ (-1);
DeltaSAeroStreamCowl(jj,1) = cp * log(TempAeroStreamCowl(jj,1)/...
Tstag) - Rgas * log(PressAeroStreamCowl(jj,1)/Pstag);
end
% Check to see if the expansion in the first expansion region is
% complete
if vAeroStream(ii) >= vRegionOne
vRegionOneCheck = 0;
else
vRegionOneCheck = 1;
end
end
% Now the program will calculate the Aerospike wall contour for the
% external portion of the nozzle. This will use C+ characteristics
% equations. The expansion point will be the last calculated cowl
% streamline point. The flow must also be turned back to Theta = 0
% Initialize loop conditions
% use ii and jj from previous loop; essentially ii = ii and jj = jj
MaContinue = 1;
DeltaVAeroTemp = 0.0;
% Loop will calculate the external aerospike wall contour
while MaContinue == 1;
ii = ii + 1; % Increase indices every iteration
%Set Expansion point conditions
ThetaAeroStreamCowl(jj,1) = ThetaAeroStreamCowl(jj,1) - DeltaVAeroTemp;
vAeroStreamCowl(jj,1) = vAeroStreamCowl(jj,1) + DeltaVAeroTemp;
% Calculate the Mach number at the current point
vRad = vAeroStreamCowl(jj,1);
MachG = MaAeroStreamCowl(jj,1);
PMtoMA % Call subroutine to find point's Mach Number
MaAeroStreamCowl(jj,1) = Mach;
AlphaAeroStreamCowl(jj,1) = asin(1/MaAeroStreamCowl(jj,1));
% Slope of the Mach line in radians
LineSlope = ThetaAeroStreamCowl(jj,1) + AlphaAeroStreamCowl(jj,1);
% Calculate the point that satisfies the Stream Function condition at
% the last expansion point on the contour
rStart = rAeroStream((ii-1),1);
% Calculate r-intercept of the Mach line
rIntercept = rAeroStreamCowl(jj,1) -...
(tan(LineSlope) * xAeroStreamCowl(jj,1));
% Initialize the values of the last point on the streamline
ThetaLast = ThetaAeroStream((ii-1),1);
xLast = xAeroStream((ii-1),1);
rLast = rAeroStream((ii-1),1);
% Calculate the point that satisfies the Stream Function
a = -tan(ThetaLast);
b = -tan(LineSlope);
c = rLast - tan(ThetaLast) * xLast;
d = rIntercept;
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rAeroStream(ii,1) = solution(1,1);
xAeroStream(ii,1) = solution(2,1);
ThetaAeroStream(ii,1) = ThetaAeroStreamCowl(jj,1);
vAeroStream(ii,1) = vAeroStreamCowl(jj,1);
MaAeroStream(ii,1) = MaAeroStreamCowl(jj,1);
AlphaAeroStream(ii,1) = AlphaAeroStreamCowl(jj,1);
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressAeroStream(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaAeroStream(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempAeroStream(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaAeroStream(ii,1)) ^ 2))) ^ (-1);
DeltaSAeroStream(ii,1) = cp * log(TempAeroStream(ii,1)/Tstag)...
- Rgas * log(PressAeroStream(ii,1)/Pstag);
end
DeltaVAeroTemp = DeltaVAeroD;
if vAeroStreamCowl(jj,1) >= vmax
MaContinue = 0;
else
MaContinue = 1;
end
end
% Program will calculate the truncated contour if the user has decided to
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
% The program will now "flip" the contour so that the axisymmetric line
% aligns with r=0.
rExit = rAeroStream(ii,1);
for ll = 1:1:(ii)
rAeroStream(ll,1) = rExit - rAeroStream(ll,1);
end
for ll = 1:1:(jj)
rAeroStreamCowl(ll,1) = rExit - rAeroStreamCowl(ll,1);
end
% "Flip" the center of the circle defining the expansion contour of Region
% One
rCenter = rExit - rCenter;
% This section "flips" the truncated Aerospike
if Truncate == 1
for jj=1:1:(zz-1)
rAeroStreamTrunc(jj,1) = rExit - rAeroStreamTrunc(jj,1); 
end
end

