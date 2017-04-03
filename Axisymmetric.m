%**************************************************************************
%* SAD Axially Symmetric Nozzle
%* Brandon Denton
%* RIT Graduate Student
%* Oct. 6, 2006
%**************************************************************************
% Program will use a combination of Shapiro, Anderson and Denton Methods to
% calculate the contour of an Axially Symmetric Supersonic Nozzle.
%**************************************************************************
%**************************************************************************
%* Directory of Variables
%**************************************************************************
%* 1. a, b, c, d ----------- Variables that store values for the matrices
%* used to calculate the coordinates and flow
%* properties of the point where the C- and C+
%* Characteristics intersect
%* 2. A, B ----------------- Matrices used the calculate the coordinates
%* and flow properties of the point where the C-
%* and C+ Characteristics intersect
%* 3. AlphaSAD() ----------- Array that holds the Mach angle of the points
%* 4. AlphaSADtemp --------- Variable that holds the temporary Mach angle
%* of the point from which the C- Characteristic
%* emanates from when the point to be calculated
%* is on the axisymmetric line
%* 5. AlphaSADtempLeft ----- Variable that temporarily stores teh Mach
%* angle of the point from which the C+
%* Characteristic emanates from
%* 6. AlphaSADtempRight ---- Variable that temporarily stores the Mach
%* angle of the point from which the C-
%* Characteristic emanates from
%* 7. AlphaSADtempRight ---- Variable that holds the temporary value of
%* the Mach angle from which the C-
%* Characteristic emanates from
%* 8. BackContinue --------- Control variable for the loop that makes sure
%* that the last point calculated on the
%* "backward" characteristics will be greater
%* than the point which will satisfy the Stream
%* Function
%* 9. Beta ----------------- Variable that holds the multipling factor of
%* the throat and dictates the radius of the
%* circle which defines the arc of the expansion
%* region
%* 10. CalcContinue --------- Control variable for the loop which
%* calculates the points that satisfy the Stream
%* Function and define the nozzle's contour
%* 11. Char ----------------- Variable which holds the number of
%* characteristics crossed when calculating the
%* nozzle's contour by the typical 2D method
%* 12. CharCheck ------------ Variable which keep tracks of how many C+
%* Characteristics are crossed during the
%* calculation of the "backward" characteristics
%* 13. DeltaR --------------- The difference between the solution to the
%* Stream Function and the value on the line
%* segment where the point lies
%* 14. DeltaTheta ----------- Control variable and flow direction error for
%* the convergence of the points position and
%* properties
%* 15. DeltaVAeroD ---------- Variable that holds the change in angle of the
%* arc that defines the expansion region
%* 16. DeltaX --------------- Variable that holds the change in x value that
%* will be used during the "backward" C-
%* Characteristic calculations
%* 17. fTheta --------------- The theta-intercept of the line segment on
%* which a point that satisfies the Stream
%* Function can be found
%* 18. fr ------------------- The r-intercept of line segment on which a
%* point which satisfies the Stream Function can
%* be found
%* 19. Gamma ---------------- Ratio of specific heats of the working fluid
%* 20. ii ------------------- Indices index
%* 21. jj ------------------- Indices index
%* 22. kk ------------------- Indices index
%* 23. ll ------------------- Indices index
%* 24. Mach ----------------- Solution value returned by the subroutine
%* PMtoMA for the Mach Number at a point
%* 25. MachG ---------------- Initial guess value used by the subroutine
%* PMtoMA to find the point's Mach Number
%* 26. MaContinue ----------- Control variable for the loop that calculates
%* the flowfield until the desired exit Mach
%* Number is reached on the axisymmetric line
%* 27. MaSAD() -------------- Array that holds the Mach Numbers of the
%* points
%* 28. MaSADtemp ------------ Variable that holds the temporary Mach Number
%* of the point from which the C- Characterisitc
%* emanates from when the point to be calculated
%* is located on the axisymmetric line
%* 29. MaSADtempLeft -------- Variable that temporarily stores the Mach
%* Number of the point from which the C+
%* Characteristic emanates from
%* 30. MaSADtempRight ------- Variable that temporarily stores the Mach
%* Number of the point from which the C-
%* Characteristic emanates from
%* 31. Mexit ---------------- User-defined desired exit Mach Number of the
%* nozzle
%* 32. NumChar -------------- Number of characteristics used in the
%* calculation
%* 33. PMtoMA --------------- Subroutine used to find the Mach Number of a
%* point in the flowfield
%* 34. PointNum ------------- Number of points in the flowfield until the
%* desired exit Mach number is reached on the
%* axisymmetric line
%* 35. Position() ----------- Array that holds the value which dictates the
%* type of flowfield point (1 = wall, 2 = flow,
%* 3 = axisymmetric)
%* 36. rCalc ---------------- The Stream Function solution of the last
%* point on the "backward" characteristic to
%* make sure it is greater than the point which
%* satisfies the Stream Function
%* 37. rCheck --------------- R-component of the last point on the
%* "backward" characteristic to check if the
%* point is greater that the point which
%* satisfies the Stream Function
%* 38. rEnd ----------------- R-component of the point which defines the
%* end of the line segment on which a point
%* which satisfies the Stream Function can be
%* found
%* 39. rLast ---------------- R-component of the point which last satisfed
%* the Stream Function
%* 40. rSAD() --------------- Array that holds the r-component of the point
%* 41. rSADcontour ---------- Variable that holds the r-component of the
%* points which define the nozzles contour
%* calculated by the 2D method
%* 42. rSADtemp ------------- Variable that holds the temporary r-component
%* of the point from which the C- Characteristic
%* emanates from when the point to be calculated
%* is located on the axisymmetric line
%* 43. rSADtempLeft --------- Variable that temporarily stores the
%* r-component of the point from which the C+
%* Characteristic is emanating from
%* 44. rSADtempRight -------- Variable that temporarily stores the
%* r-component from the point where the C-
%* Characteristic is emanating from
%* 45. rSlope --------------- The change in r for a change in x on the line
%* segment for which a point which satisfies the
%* Stream Function can be found
%* 46. rStart --------------- R-component of the point which defines the
%* start of the line segment on which the point
%* which satisfies the Stream Function can be
%* found
%* 47. rStreamContour() ----- Array which holds the r-component of the
%* points which satisify the Stream Function and
%* define the nozzle's contour
%* 48. rThroat -------------- User-defined throat radius
%* 49. solution() ----------- Solution to the Matrices used to calculate
%* the coordinates and flow properties of the
%* intersection of the C- and C+ Characteristics
%* 50. StepSize ------------- Step multiplier used when incrementally
%* increaseing the number of "backward"
%* characteristics used in the calculation
%* 51. StepSizePart --------- The decimal part of variable StepSize aslo
%* used to add in the difference between 1 and
%* the decimal part of StepSize back in to
%* StepSize to make it an integer
%* 52. StepSizeWhole -------- The whole number part of variable StepSize
%* 53. StreamContinue ------- Control variable for the loop that calculates
%* the nozzle's contour using the streamline
%* method
%* 54. ThetaCheck ----------- Flow direction of the last point on the
%* "backward" Characteristics to check if the
%* point is greater than the point which
%* satisfies the Stream Function
%* 55. ThetaEnd ------------- Flow direction of the point which defines the
%* end of the line segment which a point that
%* satisfies the Stream Function can be found
%* 56. ThetaLast ------------ Temporary value of the flow direction of the
%* last point that satisfies the Stream Function
%* also dubbes as a temporary value in the
%* convergence of DeltaTheta
%* 57. ThetaSAD() ----------- Array that holds the flow direction at the
%* point
%* 58. ThetaSADtemp --------- Variable that holds the temporary flow
%* direction of the point from which the C-
%* Characteristic emanates from when the point
%* to be calculated is located on the
%* axisymmetric line
%* 59. ThetaSADtempLeft ----- Variable that temporarily stores the flow
%* direction of the point from which the C+
%* Characteristic is emananting from
%* 60. ThetaSADtempRight ---- Variable that temporarily stores the flow
%* direction of the point where the C-
%* Characteristic is emanating from
%* 61. ThetaSlope ----------- The change in flow direction per change in x
%* of the line segment on which a point that
%* satisfies the Stream Function can be found
%* 62. ThetaStart ----------- Flow direction of the point which defines the
%* start of the line segment on which a point
%* which satisfies the Stream Function can be
%* found
%* 63. ThetaStreamContour() - Array which holds the flow direction of the
%* points which satisfy the Stream Function and
%* define the nozzle's contour
%* 64. UsedChar ------------- Temporarily holds the number of C+
%* Characteristics crossed when calculating the
%* point which satisfies the Stream Function
%* 65. vRad ----------------- Variable used by subroutine PMtoMA to find
%* the Mach Number of the point
%* 66. vSAD() --------------- Array that holds the Prandtl-Meyer expansion
%* angle of the point
%* 67. vSADtemp ------------- Variable that holds the temporary
%* Prandtl-Meyer expansion angle of the point
%* from which the C- Characteristic emanates
%* from when the point ot be calculated is
%* located on the axisymmetric line
%* 68. vSADtempLeft --------- Variable that temporarily stores the
%* Prandtl-Meyer expansion angle of the point
%* from which the C+ Characteristic emanates
%* 69. vSADtempRight -------- Variable that temporarily stores the
%* Prandtl-Meyer expansion angle of the point
%* from which the C- Characteristic emanates
%* 70. xCheck --------------- X-component of the last point on the
%* "backward" characteristics to check if it is
%* greater than the point which satisfies the
%* Stream Function
%* 71. xEnd ----------------- X-component of the point which defines the end
%* of the line segment on which a point which
%* satisfies the Stream Function can be found
%* 72. xLast ---------------- X-component of the point which last satisfied
%* the Stream Function
%* 73. xSAD() --------------- Array that holds the x-component of the point
%* 74. xSADcontour ---------- Array that holds the x-component of the
%* points which define the nozzle's contour
%* calculate by the method for 2D nozzles
%* 75. xSADtemp ------------- Variable that holds the temporary x-component
%* of the point from which the C- Characteristic
%* emanates from when the point to be calculated
%* is located on the axisymmetric line
%* 76. xSADtempLeft --------- Variable that temporarily stores the
%* x-component of the point where the C+
%* Characteristic is emanatine from
%* 77. xSADtempRight -------- Variable that temporarily stores the
%* x-component of the point where the C-
%* Characteristic is emanating from
%* 78. xStart --------------- X-component of the starting point of the line
%* segment on which the point which satisfies
%* the Stream Function can be found
%* 79. xStreamContour() ----- Array which holds the x-components of the
%* points which satisfy the Stream Function and
%* define the nozzle's contour
%* 80. zz ------------------- Indices index
%**************************************************************************
%**************************************************************************
%* Start Program
%**************************************************************************
%**************************************************************************
format long
DeltaX = DeltaVAeroD;
%The calculation must be started with the first characteristic
NumChar = 1;
PointNum = 0; %Initializes the variable PointNum
ii = 1; %Initializes the looping variable ii
jj = 1; %Initializes the looping variable jj
MaContinue = 0;
while MaContinue == 0
PointNum = PointNum + (NumChar + 1);
%Initialize the position type of the points
for ii = (PointNum - NumChar):1:PointNum
if ii == (PointNum - NumChar)
Position(ii,1) = 1;
elseif ii == PointNum
Position(ii,1) = 3;
else
Position(ii,1) = 2;
end
end
for ii = (PointNum - NumChar):1:PointNum
if Position(ii,1) == 1
%Set the angle at which the calculation has been swept through
ThetaSAD(ii,1) = NumChar * DeltaVAeroD;
%Calculate the x-coordinate of the point that corresponds to the
%current characteristic (Wall point)
xSAD(ii,1) = (Beta * rThroat) * sin(ThetaSAD(ii,1));
%Calculate the radial position of the point on the wall of the
%current characteristic (Wall point)
rSAD(ii,1) = rThroat + ((Beta * rThroat) *...
(1 - cos(ThetaSAD(ii,1))));
%Calculate the Prandtl-Meyer Expansion Angle vSAD() for the
%current point (Wall point)
vSAD(ii,1) = ThetaSAD(ii,1);
%Now calculate the Mach number and Alpha associated with the
%point. Must convert for subprogram to find the associated Mach
%Number
vRad = vSAD(ii,1);
MachG = 1.0;
PMtoMA %calls subprogram to find Mach Number
MaSAD(ii,1) = Mach;
AlphaSAD(ii,1) = asin(1/MaSAD(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point log = natural log
PressSAD(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempSAD(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-1);
DeltaS(ii,1) = cp * log(TempSAD(ii,1)/Tstag) - Rgas *...
log(PressSAD(ii,1)/Pstag);
end
elseif Position(ii,1) == 3
%Calculate the position of the point when it is on the
%axisymmetric line
ThetaSAD(ii,1) = 0;
rSAD(ii,1) = 0;
%Initialize temperary values for the loop calculation of the
%curved characteristics
ThetaSADtemp = ThetaSAD((ii-1),1);
rSADtemp = rSAD((ii-1),1);
xSADtemp = xSAD((ii-1),1);
AlphaSADtemp = AlphaSAD((ii-1),1);
MaSADtemp = MaSAD((ii-1),1);
vSADtemp = vSAD((ii-1),1);
DeltaTheta = 1.0; %Initialize the loop below
ThetaLast = 100;
while DeltaTheta >= 1e-10
%Calculate the position of the current point
xSAD(ii,1) = ((rSADtemp - (tan(ThetaSADtemp -...
AlphaSADtemp) * xSADtemp)) /...
(-tan(ThetaSADtemp - AlphaSADtemp))) -...
rSAD(ii,1);
%Calculate the flow properties of the current point
vSAD(ii,1) = ThetaSADtemp + vSADtemp - ThetaSAD(ii,1) + ...
((1/(sqrt((MaSADtemp^2) - 1) -...
(1/tan(ThetaSADtemp)))) * ((rSAD(ii,1) -...
rSADtemp) / rSADtemp));
%Calculate the Mach Number at the current point
vRad = vSAD(ii,1);
MachG = 1.0;
PMtoMA %calls subprogram to find Mach Number
MaSAD(ii,1) = Mach;
AlphaSAD(ii,1) = asin(1/MaSAD(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point
PressSAD(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempSAD(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-1);
DeltaS(ii,1) = cp * log(TempSAD(ii,1)/Tstag) - Rgas*...
log(PressSAD(ii,1)/Pstag);
end
%Calculate the change in Theta for the loop
DeltaTheta = abs(ThetaLast - ThetaSAD(ii,1));
if DeltaTheta > 1e-10
ThetaLast = ThetaSAD(ii,1);
% Calculate averages of all values and replace the
% "temp" values with these. This is an approximation
% that the charactertistics are curved
ThetaSADtemp = (ThetaSADtemp + ThetaSAD(ii,1)) / 2;
AlphaSADtemp = (AlphaSADtemp + AlphaSAD(ii,1)) / 2;
MaSADtemp = (MaSADtemp + MaSAD(ii,1)) / 2;
rSADtemp = (rSADtemp + rSAD(ii,1)) / 2;
xSADtemp = (xSADtemp + xSAD(ii,1)) / 2;
vSADtemp = (vSADtemp + vSAD(ii,1)) / 2;
end
end
%Check to see if the point on the axisymmetric point has
%achieved the desired exit mach number
if MaSAD(ii,1) >= Mexit
MaContinue = 1;
else
MaContinue = 0;
end
%Increase the Number of Characteristics
NumChar = NumChar + 1;
else
%Initialize C- Characteristic values for the calculation of the
%flowfield point
ThetaSADtempRight = ThetaSAD((ii-1),1);
rSADtempRight = rSAD((ii-1),1);
xSADtempRight = xSAD((ii-1),1);
AlphaSADtempRight = AlphaSAD((ii-1),1);
MaSADtempRight = MaSAD((ii-1),1);
vSADtempRight = vSAD((ii-1),1);
%Initialize C+ Characteristic values for the calculation of the
%flowfield point
ThetaSADtempLeft = ThetaSAD((ii-NumChar),1);
rSADtempLeft = rSAD((ii-NumChar),1);
xSADtempLeft = xSAD((ii-NumChar),1);
AlphaSADtempLeft = AlphaSAD((ii-NumChar),1);
MaSADtempLeft = MaSAD((ii-NumChar),1);
vSADtempLeft = vSAD((ii-NumChar),1);
DeltaTheta = 1.0; %Initialize the following loop
ThetaLast = 100;
while DeltaTheta >= 1e-10
%Calculate the position of the current point
a = tan(ThetaSADtempRight - AlphaSADtempRight);
b = tan(ThetaSADtempLeft + AlphaSADtempLeft);
c = rSADtempRight - (a * xSADtempRight);
d = rSADtempLeft - (b * xSADtempLeft);
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
rSAD(ii,1) = solution(1,1);
xSAD(ii,1) = solution(2,1);
%Calculate the flow properties of the current point
c = (ThetaSADtempRight + vSADtempRight) + ((1 /...
((sqrt(MaSADtempRight^2 - 1)) -...
(1/tan(ThetaSADtempRight)))) * ((rSAD(ii,1) -...
rSADtempRight) / rSADtempRight));
if ThetaSADtempLeft == 0.0
d = (2 * ThetaSADtempLeft) - vSADtempLeft;
A = [1 1; 2 -1];
else
d = (ThetaSADtempLeft - vSADtempLeft) - ((1 /...
((sqrt(MaSADtempLeft^2 - 1)) +...
(1/tan(ThetaSADtempLeft)))) * ((rSAD(ii,1) -...
rSADtempLeft) / rSADtempLeft));
A = [1 1; 1 -1];
end
B = [c; d];
solution = A\B;
ThetaSAD(ii,1) = solution(1,1);
vSAD(ii,1) = solution(2,1);
%Calculate the Mach Number at the current point
vRad = vSAD(ii,1);
MachG = 1.0;
PMtoMA %calls subprogram to find the Mach Number
MaSAD(ii,1) = Mach;
AlphaSAD(ii,1) = asin(1/MaSAD(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point
PressSAD(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempSAD(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-1);
DeltaS(ii,1) = cp * log(TempSAD(ii,1)/Tstag) - Rgas*...
log(PressSAD(ii,1)/Pstag);
end
%Calculate the change in Theta for the loop
DeltaTheta = abs(ThetaLast - ThetaSAD(ii,1));
if DeltaTheta > 1e-10
ThetaLast = ThetaSAD(ii,1);
%Calculate averages of all values and replace the
%"temp" values with these. This is an approximation
%that the charactertistics are curved.
%C- Characteristic
ThetaSADtempRight = (ThetaSADtempRight +...
ThetaSAD(ii,1)) / 2;
AlphaSADtempRight = (AlphaSADtempRight +...
AlphaSAD(ii,1)) / 2;
MaSADtempRight = (MaSADtempRight + MaSAD(ii,1)) / 2;
rSADtempRight = (rSADtempRight + rSAD(ii,1)) / 2;
xSADtempRight = (xSADtempRight + xSAD(ii,1)) / 2;
vSADtempRight = (vSADtempRight + vSAD(ii,1)) / 2;
%C+ Characteristic
ThetaSADtempLeft = (ThetaSADtempLeft +...
ThetaSAD(ii,1)) / 2;
AlphaSADtempLeft = (AlphaSADtempLeft +...
AlphaSAD(ii,1)) / 2;
MaSADtempLeft = (MaSADtempLeft + MaSAD(ii,1)) / 2;
rSADtempLeft = (rSADtempLeft + rSAD(ii,1)) / 2;
xSADtempLeft = (xSADtempLeft + xSAD(ii,1)) / 2;
vSADtempLeft = (vSADtempLeft + vSAD(ii,1)) / 2;
end
end
end
end
end
%**************************************************************************
%**************************************************************************
%Need to calculate the the points on the contour of the nozzle. We will
%assume that at straight line eminates from the last point on the entrance
%region with a slope of the average of the point and the next point on the
%last calculated characteristic. We will assume that the fluid at the wall
%will have the same properties as the point on the last calculated
%characteristic and the C+ characteristic from this point and the previous
%wall point intersect will be the wall contour point. (This may sound
%confusing)
%Find the first wall point on the wall contour
ii = PointNum + 1; %Initialize the ii loop index and the index of the first
%point on the wall contour
%Calculate the first wall contour point
ThetaSAD(ii,1) = ThetaSAD((ii-(NumChar-1)),1);
a = tan(ThetaSAD((ii-(NumChar-1)),1) + AlphaSAD((ii-(NumChar-1)),1));
b = (ThetaSAD((ii-NumChar),1) + ThetaSAD((ii-(NumChar-1)),1)) / 2;
c = rSAD((ii-(NumChar-1)),1) - (a * xSAD((ii-(NumChar-1)),1));
d = rSAD((ii-NumChar),1) - (b * xSAD((ii-NumChar),1));
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
rSAD(ii,1) = solution(1,1);
xSAD(ii,1) = solution(2,1);
%Calculate the rest of the wall contour points
for jj = (ii + 1):1:(PointNum+(NumChar-1))
ThetaSAD(jj,1) = ThetaSAD((jj-(NumChar-1)),1);
a = tan(ThetaSAD((jj-(NumChar-1)),1) + AlphaSAD((jj-(NumChar-1)),1));
b = (ThetaSAD((jj-1),1) + ThetaSAD((jj-(NumChar-1)),1)) / 2;
c = rSAD((jj-(NumChar-1)),1) - (a * xSAD((jj-(NumChar-1)),1));
d = rSAD((jj-1),1) - (b * xSAD((jj-1),1));
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
rSAD(jj,1) = solution(1,1);
xSAD(jj,1) = solution(2,1);
end
%**************************************************************************
%**************************************************************************
%Program will now retain only the points that lie on the contour of the
%wall
jj = 1; %Initialize the index of finding the right point on the wall
Char = 1; %Initialize the variable that holds the temperary characteristic
%number
for ii = 1:1:(NumChar-1)
xSADcontour(ii,1) = xSAD(jj,1);
rSADcontour(ii,1) = rSAD(jj,1);
jj = jj + (Char + 1);
Char = Char + 1;
end
jj = PointNum + 1;
for ii = NumChar:1:((NumChar-1)+(NumChar-1))
xSADcontour(ii,1) = xSAD(jj,1);
rSADcontour(ii,1) = rSAD(jj,1);
jj = jj + 1;
end
%**************************************************************************
%**************************************************************************
%Another approach to calculate for the contour that uses the streamline and
%"backward" calculated characteristics
StreamContinue = 1;
%Initialize Stepsize to reduce the number of wasted iterations of the prog.
%StepSize = PointNum/(NumChar-1);
%StepSizeWhole = fix(StepSize);
%StepSizePart = StepSize - StepSizeWhole;
%StepSizePart = 1 - StepSizePart;
%StepSize = StepSize + StepSizePart + 1;
StepSize = 1.0;
while StreamContinue == 1
kk = PointNum + 1;
ThetaSAD(kk,1) = ThetaSAD((kk-1),1);
vSAD(kk,1) = vSAD((kk-1),1);
xSAD(kk,1) = xSAD((kk-1),1) + (StepSize * DeltaX);
MaSAD(kk,1) = MaSAD((kk-1),1);
AlphaSAD(kk,1) = AlphaSAD((kk-1),1);
%Calculate the r position of the next point
rSAD(kk,1) = rSAD((kk-1),1) + (tan(ThetaSAD((kk-1),1) +...
AlphaSAD((kk-1),1)) * (xSAD(kk,1) - xSAD((kk-1),1)));
%Calculate all of the points along the first "backward" characteristic
zz = 3;
for ii = (PointNum+2):1:(PointNum+NumChar)
%Initialize C- Characteristic values for the calculation of the
%flowfield point
ThetaSADtempRight = ThetaSAD((ii-1),1);
rSADtempRight = rSAD((ii-1),1);
xSADtempRight = xSAD((ii-1),1);
AlphaSADtempRight = AlphaSAD((ii-1),1);
MaSADtempRight = MaSAD((ii-1),1);
vSADtempRight = vSAD((ii-1),1);
%Initialize C+ Characteristic values for the calculation of the
%flowfield point
ThetaSADtempLeft = ThetaSAD((ii-zz),1);
rSADtempLeft = rSAD((ii-zz),1);
xSADtempLeft = xSAD((ii-zz),1);
AlphaSADtempLeft = AlphaSAD((ii-zz),1);
MaSADtempLeft = MaSAD((ii-zz),1);
vSADtempLeft = vSAD((ii-zz),1);
DeltaTheta = 1.0; %Initialize the following loop
ThetaLast = 100;
while DeltaTheta >= 1e-10
%Calculate the position of the current point
a = tan(ThetaSADtempRight - AlphaSADtempRight);
b = tan(ThetaSADtempLeft + AlphaSADtempLeft);
c = rSADtempRight - (a * xSADtempRight);
d = rSADtempLeft - (b * xSADtempLeft);
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
rSAD(ii,1) = solution(1,1);
xSAD(ii,1) = solution(2,1);
%Calculate the flow properties of the current point
if ThetaSADtempRight == 0
c = ThetaSADtempRight + vSADtempRight;
A = [1 1; 1 -1];
else
c = (ThetaSADtempRight + vSADtempRight) + ((1 /...
((sqrt(MaSADtempRight^2 - 1)) -...
(1/tan(ThetaSADtempRight)))) * ((rSAD(ii,1) -...
rSADtempRight) / rSADtempRight));
end
d = (ThetaSADtempLeft - vSADtempLeft) - ((1 /...
((sqrt(MaSADtempLeft^2 - 1)) +...
(1/tan(ThetaSADtempLeft)))) * ((rSAD(ii,1) -...
rSADtempLeft) / rSADtempLeft));
A = [1 1; 1 -1];
B = [c; d];
solution = A\B;
ThetaSAD(ii,1) = solution(1,1);
vSAD(ii,1) = solution(2,1);
%Calculate the Mach Number at the current point
vRad = vSAD(ii,1);
MachG = 1.0;
PMtoMA %calls subprogram to find the Mach Number
MaSAD(ii,1) = Mach;
AlphaSAD(ii,1) = asin(1/MaSAD(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point
PressSAD(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempSAD(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-1);
DeltaS(ii,1) = cp * log(TempSAD(ii,1)/Tstag) - Rgas *...
log(PressSAD(ii,1)/Pstag);
end
%Calculate the change in Theta for the loop
DeltaTheta = abs(ThetaLast - ThetaSAD(ii,1));
if DeltaTheta > 1e-10
ThetaLast = ThetaSAD(ii,1);
%Calculate averages of all values and replace the "temp"
%values with these. This is an approximation that the
%charactertistics are curved.
%C- Characteristic
ThetaSADtempRight = (ThetaSADtempRight +...
ThetaSAD(ii,1)) / 2;
AlphaSADtempRight = (AlphaSADtempRight +...
AlphaSAD(ii,1)) / 2;
MaSADtempRight = (MaSADtempRight + MaSAD(ii,1)) / 2;
rSADtempRight = (rSADtempRight + rSAD(ii,1)) / 2;
xSADtempRight = (xSADtempRight + xSAD(ii,1)) / 2;
vSADtempRight = (vSADtempRight + vSAD(ii,1)) / 2;
%C+ Characteristic
ThetaSADtempLeft = (ThetaSADtempLeft + ThetaSAD(ii,1)) / 2;
AlphaSADtempLeft = (AlphaSADtempLeft + AlphaSAD(ii,1)) / 2;
MaSADtempLeft = (MaSADtempLeft + MaSAD(ii,1)) / 2;
rSADtempLeft = (rSADtempLeft + rSAD(ii,1)) / 2;
xSADtempLeft = (xSADtempLeft + xSAD(ii,1)) / 2;
vSADtempLeft = (vSADtempLeft + vSAD(ii,1)) / 2;
end
end
zz = zz + 2;
end
% This if statement makes sure that the last point calculated in the
% flow has an r-coordinate greater than the r-coordinate of the last
% expansion point. If it doesn't the stream function solution will
% fail.
if rSAD((PointNum+NumChar),1) > rSAD((PointNum-(NumChar-1)),1)
StreamContinue = 0;
else
StreamContinue = 1;
StepSize = StepSize + 1;
end
end
%Calculate the point that satisfies the Stream Function condition at the
%last expansion point on the contour
xStart = xSAD((ii-1),1);
rStart = rSAD((ii-1),1);
ThetaStart = ThetaSAD((ii-1),1);
ThetaEnd = ThetaSAD(ii,1);
rEnd = rSAD(ii,1);
xEnd = xSAD(ii,1);
%Calculate the slope and y-intercept of the straight line that is
%approximating the characteristic between the last two calculated
%characteristic points
A = [xStart 1; xEnd 1];
B = [rStart; rEnd];
solution = A\B;
rSlope = solution(1,1); %Slope of r with respect to x
fr = solution(2,1); %y-intercept of the r line
B = [ThetaStart; ThetaEnd];
solution = A\B;
ThetaSlope = solution(1,1); %Slope of Theta with respect to x
fTheta = solution(2,1); %y-intercept of the Theta line
%Initiate values of the last point on the streamline
ThetaLast = ThetaSAD((ii-((2*NumChar)-1)),1);
xLast = xSAD((ii-((2*NumChar)-1)),1);
rLast = rSAD((ii-((2*NumChar)-1)),1);
%Calculate the x- and Theta- coordinate on the line based on the
%r-coordinate of the previous point that satisfies the stream function.
%This is avoid extra computational time due to the fact that the point next
%point that satisfies the stream function must have an r-coordinate greater
%than the previous stream line point's r-coordinate
%Calculate the position of the point the satisfies the streamline condition
a = -tan(ThetaLast);
b = -rSlope;
c = rLast - tan(ThetaLast) * xLast;
d = fr;
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rStreamContour(1,1) = solution(1,1);
xStreamContour(1,1) = solution(2,1);
ThetaStreamContour(1,1) = ThetaSlope * xStreamContour(1,1) + fTheta;
%Calculate the rest of the streamline points by continuing the calculation
%until the number of characteristics crossed is equal to 1
ll = 2; %Index for stream line points
UsedChar = NumChar;
CalcContinue = 1;
if MaSAD((PointNum),1) > 4.5
StepSize = (2 * StepSize);
else
Stepsize = StepSize;
end
while CalcContinue ==1
CharCheck = 1; %Resets the number of Characteristics calculated
kk = ii + 1;
ThetaSAD(kk,1) = ThetaSAD((kk-UsedChar),1);
vSAD(kk,1) = vSAD((kk-UsedChar),1);
xSAD(kk,1) = xSAD((kk-UsedChar),1) + ((StepSize) * DeltaX);
MaSAD(kk,1) = MaSAD((kk-UsedChar),1);
AlphaSAD(kk,1) = AlphaSAD((kk-UsedChar),1);
%Calculate the r position of the next point
rSAD(kk,1) = rSAD((kk-UsedChar),1) + (tan(ThetaSAD((kk-UsedChar),1)+...
AlphaSAD((kk-UsedChar),1)) * (xSAD(kk,1) - xSAD((kk-UsedChar),1)));
if rSAD(kk,1) > rStreamContour((ll-1),1)
CalcContinue = 0;
else
ii = ii + 1; %Increases the ii index to match kk index
%Calculate all of the points along the first "backward" characteristic
BackContinue = 1.0;
while BackContinue == 1.0
ii = ii + 1; %Increases the point index
%Initialize C- Characteristic values for the calculation of the
%flowfield point
ThetaSADtempRight = ThetaSAD((ii-1),1);
rSADtempRight = rSAD((ii-1),1);
xSADtempRight = xSAD((ii-1),1);
AlphaSADtempRight = AlphaSAD((ii-1),1);
MaSADtempRight = MaSAD((ii-1),1);
vSADtempRight = vSAD((ii-1),1);
%Initialize C+ Characteristic values for the calculation of the
%flowfield point
ThetaSADtempLeft = ThetaSAD((ii-UsedChar),1);
rSADtempLeft = rSAD((ii-UsedChar),1);
xSADtempLeft = xSAD((ii-UsedChar),1);
AlphaSADtempLeft = AlphaSAD((ii-UsedChar),1);
MaSADtempLeft = MaSAD((ii-UsedChar),1);
vSADtempLeft = vSAD((ii-UsedChar),1);
DeltaTheta = 1.0; %Initialize the following loop
ThetaLast = 100;
while DeltaTheta >= 1e-10
%Calculate the position of the current point
a = tan(ThetaSADtempRight - AlphaSADtempRight);
b = tan(ThetaSADtempLeft + AlphaSADtempLeft);
c = rSADtempRight - (a * xSADtempRight);
d = rSADtempLeft - (b * xSADtempLeft);
A = [1 -a; 1 -b];
B = [c; d];
solution = A\B;
rSAD(ii,1) = solution(1,1);
xSAD(ii,1) = solution(2,1);
%Calculate the flow properties of the current point
if ThetaSADtempRight == 0
c = ThetaSADtempRight + vSADtempRight;
else
c = (ThetaSADtempRight + vSADtempRight) + ((1 /...
((sqrt(MaSADtempRight^2 - 1)) -...
(1/tan(ThetaSADtempRight)))) * ((rSAD(ii,1) -...
rSADtempRight) / rSADtempRight));
end
d = (ThetaSADtempLeft - vSADtempLeft) - ((1 /...
((sqrt(MaSADtempLeft^2 - 1)) +...
(1/tan(ThetaSADtempLeft)))) * ((rSAD(ii,1) -...
rSADtempLeft) / rSADtempLeft));
A = [1 1; 1 -1];
B = [c; d];
solution = A\B;
ThetaSAD(ii,1) = solution(1,1);
vSAD(ii,1) = solution(2,1);
%Calculate the Mach Number at the current point
vRad = vSAD(ii,1);
MachG = 1.0;
PMtoMA %calls subprogram to find the Mach Number
MaSAD(ii,1) = Mach;
AlphaSAD(ii,1) = asin(1/MaSAD(ii,1));
if EntroCheck == 1
% Calculate the Pressure, Temperature and change in Entropy
% at the point
PressSAD(ii,1) = Pstag * ((1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-Gamma/(Gamma - 1)));
TempSAD(ii,1) = Tstag * (1 + (((Gamma - 1)/2) *...
((MaSAD(ii,1)) ^ 2))) ^ (-1);
DeltaS(ii,1) = cp * log(TempSAD(ii,1)/Tstag) - Rgas*...
log(PressSAD(ii,1)/Pstag);
end
%Calculate the change in Theta for the loop
DeltaTheta = abs(ThetaLast - ThetaSAD(ii,1));
if DeltaTheta > 1e-10
ThetaLast = ThetaSAD(ii,1);
%Calculate averages of all values and replace the
%"temp" values with these. This is an approximation
%that the charactertistics are curved.
%C- Characteristic
ThetaSADtempRight = (ThetaSADtempRight +...
ThetaSAD(ii,1)) / 2;
AlphaSADtempRight = (AlphaSADtempRight +...
AlphaSAD(ii,1)) / 2;
MaSADtempRight = (MaSADtempRight + MaSAD(ii,1)) / 2;
rSADtempRight = (rSADtempRight + rSAD(ii,1)) / 2;
xSADtempRight = (xSADtempRight + xSAD(ii,1)) / 2;
vSADtempRight = (vSADtempRight + vSAD(ii,1)) / 2;
%C+ Characteristic
ThetaSADtempLeft = (ThetaSADtempLeft +...
ThetaSAD(ii,1)) / 2;
AlphaSADtempLeft = (AlphaSADtempLeft +...
AlphaSAD(ii,1)) / 2;
MaSADtempLeft = (MaSADtempLeft + MaSAD(ii,1)) / 2;
rSADtempLeft = (rSADtempLeft + rSAD(ii,1)) / 2;
xSADtempLeft = (xSADtempLeft + xSAD(ii,1)) / 2;
vSADtempLeft = (vSADtempLeft + vSAD(ii,1)) / 2;
end
end
xLast = xStreamContour((ll-1),1);
rLast = rStreamContour((ll-1),1);
ThetaLast = ThetaStreamContour((ll-1),1);
rCheck = rSAD(ii,1);
xCheck = xSAD(ii,1);
ThetaCheck = ThetaSAD(ii,1);
rCalc = tan(ThetaLast) * (xCheck - xLast) + rLast;
DeltaR = rCheck - rCalc;
if DeltaR > 0.0
BackContinue = 0.0;
else
BackContinue = 1.0;
end
CharCheck = CharCheck + 1; %Increases the Number of Used Chars
end
UsedChar = CharCheck; % Sets UsedChar to the correct value of the
% distance of points between any two
% "backward" calculated characteristics
%Calculate the point that satisfies the Stream Function condition
%at the last expansion point on the contour
xStart = xSAD((ii-1),1);
rStart = rSAD((ii-1),1);
ThetaStart = ThetaSAD((ii-1),1);
xEnd = xSAD(ii,1);
rEnd = rSAD(ii,1);
ThetaEnd = ThetaSAD(ii,1);
%Calculate the slope and y-intercept of the straight line that is
%approximating the characteristic between the last two calculated
%characteristic points
A = [xStart 1; xEnd 1];
B = [rStart; rEnd];
solution = A\B;
rSlope = solution(1,1); %Slope of r with respect to x
fr = solution(2,1); %y-intercept of the r line
B = [ThetaStart; ThetaEnd];
solution = A\B;
ThetaSlope = solution(1,1); %Slope of Theta with respect to x
fTheta = solution(2,1); %y-intercept of the Theta line
%Initiate values of the last point on the streamline
ThetaLast = ThetaStreamContour((ll-1),1);
xLast = xStreamContour((ll-1),1);
rLast = rStreamContour((ll-1),1);
%Calculate the x- and Theta- coordinate on the line based on the
%r-coordinate of the previous point that satisfies the stream
%function. This is avoid extra computational time due to the fact
%that the point next point that satisfies the stream function must
%have an r-coordinate greater than the previous stream line point's
%r-coordinate
%Calculate the position of the point the satisfies the streamline
%condition.
a = -tan(ThetaLast);
b = -rSlope;
c = rLast - tan(ThetaLast) * xLast;
d = fr;
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rStreamContour(ll,1) = solution(1,1);
xStreamContour(ll,1) = solution(2,1);
ThetaStreamContour(ll,1) = ThetaSlope*xStreamContour(ll,1)+fTheta;
ll = ll + 1; % Increases index variable
CalcContinue = 1;
end
end
%Calculate the point that satisfies the Stream Function condition at the
%last expansion point on the contour
xStart = xSAD((kk-UsedChar),1);
rStart = rSAD((kk-UsedChar),1);
ThetaStart = ThetaSAD((kk-UsedChar),1);
xEnd = xSAD(kk,1);
rEnd = rSAD(kk,1);
ThetaEnd = ThetaSAD(kk,1);
%Calculate the slope and y-intercept of the straight line that is
%approximating the characteristic between the last two calculated
%characteristic points
A = [xStart 1; xEnd 1];
B = [rStart; rEnd];
solution = A\B;
rSlope = solution(1,1); %Slope of r with respect to x
fr = solution(2,1); %y-intercept of the r line
B = [ThetaStart; ThetaEnd];
solution = A\B;
ThetaSlope = solution(1,1); %Slope of Theta with respect to x
fTheta = solution(2,1); %y-intercept of the Theta line
%Initiate values of the last point on the streamline
ThetaLast = ThetaStreamContour((ll-1),1);
xLast = xStreamContour((ll-1),1);
rLast = rStreamContour((ll-1),1);
%Calculate the position of the point the satisfies the streamline
%condition
a = -tan(ThetaLast);
b = -rSlope;
c = rLast - tan(ThetaLast) * xLast;
d = fr;
A = [1 a; 1 b];
B = [c; d];
solution = A\B;
rStreamContour(ll,1) = solution(1,1);
xStreamContour(ll,1) = solution(2,1);
ThetaStreamContour(ll,1) = ThetaSlope * xStreamContour(ll,1) + fTheta;