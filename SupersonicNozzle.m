%**************************************************************************
%* Supersonic Nozzle Design
%* By: Brandon Denton
%* RIT Graduate Student
%**************************************************************************
%**************************************************************************
%* Dictionary of Variable
%**************************************************************************
%* 1. AxisymmetricSAD ---------- Subroutine that calculates a contour for
%* an annular supersonic nozzle using the
%* stream function and a combination of
%* techniques outlined by Anderson, JD, and
%* Shapiro, complied by Brandon Denton
%* (Thesis Attached)
%* 2. Beta --------------------- Multiplication factor for radius of
%* entrance region
%* 3. Choice ------------------- Holds the value of the Choice variable
%* for the the truncation of the aerospike
%* nozzle's contour
%* 4. DentonAerospike ---------- Subroutine that calculates the contour of
%* a supersonic aerospike nozzle written by
%* Brandon Denton, 2007. (Attached Thesis)
%* Using the stream function
%* 5. Gamma -------------------- Ratio of Specific Heats of working fluid
%* 6. Mexit -------------------- Desired exit Mach Number
%* 7. nb ----------------------- Variable that tells the subroutine
%* AnnularAerospike what percentatge of the
%* length of the idealized Aerospike nozzle
%* should be retained
%* 8. NumChar ------------------ Number of Characteristics used in the
%* Method of Characteristic Calculations
%* aerospikes nozzle's contour
%* 9. Percent ------------------ Variable that holds the percentage the
%* user would like to retain of the
%* idealized Aerospike nozzle
%* 10. Pstag -------------------- Chamber Stagnation Pressure
%* 11. rAeroExpansion ----------- Array that holds the r-component of the
%* expansion point calculated by subroutine
%* DentonAerospike
%* 12. rAeroStream -------------- Array that holds the r-component of the
%* Internal-External Aerospike nozzle
%* contour
%* 13. rAeroStreamContour ------- Array that holds the r-component of the
%* aerospike nozzle calculated by subroutine
%* DentonAerospike
%* 14. rAeroStreamCowl ---------- Array that holds the r-components of the
%* cowl of the Internal-External Aerospike
%* nozzle
%* 15. rAeroStreamContourTrunc -- Array that holds the r-component of the
%* annular aerospike calculated by
%* subroutine DentonAerospike
%* 16. rAeroStreamTrunc --------- Array that holds the r-components of the
%* truncated Internal-External Aerospike
%* nozzle
%* 17. rSADcontour -------------- Array that holds the r-component of the
%* annular nozzle contour calculated by
%* subroutine AxisymmetricSAD but calculated
%* it using the same technique used for the
%* 2D nozzles
%* 18. rStreamContour ----------- Array that holds the r-component of the
%* annular nozzle contour calculated by
%* subroutine AxisymmetricSAD using the
%* stream function technique
%* 19. rThroat ------------------ User-defined radius of the throat
%* 20. Truncate ----------------- Variable that tells the program whether
%* or not the user would like to truncate
%* the aerospike nozzle
%* 21. Tstag -------------------- Chamber Stagnation Temperature
%* 22. xAeroExpansion ----------- Array that holds the x-component of the
%* expansion point calculated by subroutine
%* DentonAerospike
%* 23. xAeroStream -------------- Array that holds the x-component of the
%* Internal-External Aerospike nozzle
%* 24. xAeroStreamContour ------- Array that holds the x-component of the
%* aerospike nozzle calculated by subroutine
%* DentonAerospike
%* 25. xAeroStreamContourTrunc -- Array that holds the x-component of the
%* annular aerospike calculated by
%* subroutine DentonAerospike
%* 26. xAeroStreamCowl ---------- Array that holds the x-components of the
%* cowl of the Internal-External Aerospike
%* nozzle
%* 27. xAeroStreamTrunc --------- Array that holds the x-component of the
%* truncated Internal-External Aerospike
%* nozzle
%* 28. xSADcontour -------------- Array that holds the x-component of the
%* annular nozzle contour calculated by
%* subroutine AxisymmetricSAD but calculated
%* it using the same technique used for the
%* 2D nozzles
%* 29. xStreamContour ----------- Array that holds the x-component of the
%* annular nozzle contour calculated by
%* subroutine AxisymmetricSAD using the
%* stream function technique
%* 30. zz ----------------------- Indices indicator used in subroutine
%* Annular Aerospike
%**************************************************************************
%**************************************************************************
%* Start Program
%**************************************************************************
%**************************************************************************
disp(' ');
% Program will ask for the desired Exit Mach Number
Mexit = input('Please enter the desired Exit Mach Number: ');
disp(' ');
% Program will ask for the ratio of specific heats for the fluid
Gamma = input('Please enter the ratio of specific heats for the fluid: ');
disp(' ');
% Program will ask if the designer would like to perform a entropy check to
% see if the change in entropy is zero or at least very small to validate
% the isentropic assumption used in the calculation
EntroLoop = 1;
while EntroLoop == 1
disp('Would you like to perform an entropy check of the flowfield?');
disp(' 1 = Yes');
disp(' 2 = No');
EntroCheck = input('-> ');
disp(' ');
% Program asks for the Rocket Chambers Pressure and Temperature
%(Stagnation or Total) and specific heat at constant pressure
if EntroCheck == 1
Pstag = input...
('Please enter the Rocket Chamber Stagnation Pressure:(Pa) ');
disp(' ');
Tstag = input...
('Please enter the Rocket Chamber Stagnation Temperature:(K) ');
disp(' ');
cp = input...
('Please enter the specific heat at constant pressure of the fluid:(J/kg-K) ');
disp(' ');
Rgas = (cp * (Gamma - 1))/Gamma;
EntroLoop = 0;
elseif EntroCheck == 2
EntroLoop = 0;
else
disp('You have entered an invalid value');
disp(' ');
end
end
% Program will ask for the radius of the throat
disp('Please enter the radius of the throat of the nozzle.');
disp('for a dimensionless nozzle, enter a radius of 1');
rThroat = input('-> ');
disp(' ');
% Program will ask for the Multiplication factor of the throat to determine
% the radius of the circle defining the entrance region
disp('FOR THE AXISYMMETRIC AND INTERNAL-EXTERNAL AEROSPIKE NOZZLE');
disp('Please enter the factor of the throat that will define the');
Beta = input...
('radius of the circle for the entrance region of the nozzles: ');
disp(' ');
% Program will ask if the user would like to calculate the contour of a
% truncated Aerospike nozzle and if so what percent of the idealized length
% would the user like to retain
zz = 0;
kk = 0;
while zz == 0
disp('Do you want to truncate the Aerospike Nozzle?');
disp(' 1 = Yes');
disp(' 2 = No');
Truncate = input('-> ');
if Truncate == 2
nb = 0.0;
zz = 1;
elseif Truncate == 1
disp(' ');
Percent = input...
('Please enter the percentage of the idealize length you would like to retain: ');
nb = Percent/100;
zz = 1;
else
disp(' ');
disp('You have entered and invalid response, please retry.');
end
end
% Program will ask for the desired change in Prandtl-Meyer Expansion Angle
% for the Aerospike design using the streamline technique
disp(' ');
DeltaVAeroD = input...
('Please enter the desired change in Prandtl-Meyer Expansion Angle: ');
% Program will now calculate the Non-Dimensional Axi-symmetric Aerospike
% Nozzle for the desired Exit Mach Number
% NOTE******* This aerospike nozzle is calculated using the streamline
% technique
disp('Calling External Aerospike...');
DentonAerospike
% Program will now calculate the Non-Dimensional Axi-symmetric
% Internal-External Aerospike nozzle for the desired Exit Mach Number
disp('Calling Internal-External Aerospike...');
IEaerospike
% Program will now calculate the Non-Dimensional Aerospike nozzle for the
% desired Exit Mach Number using the Shapiro, Anderson, Denton Method
disp('Calling Ideal Bell...');
%Axisymmetric
%disp('Back from Ideal Bell...');
% Save output to file
%save file1.dat xSADcontour
%save file2.dat rSADcontour
%save file3.dat xStreamContour
%save file4.dat rStreamContour
%save file5.dat xAeroStreamContour
%save file6.dat rAeroStreamContour
%save file7.dat xAeroExpansion
%save file8.dat rAeroExpansion
save ramp.dat xAeroStream rAeroStream
save cowl.dat xAeroStreamCowl rAeroStreamCowl
%if Truncate == 1
%  save file13.dat xAeroStreamContourTrunc
%  save file14.dat rAeroStreamTrunc
%end

% Output Plot
%if Truncate == 1
%plot(xSADcontour, rSADcontour, 'g.', xStreamContour,...
%rStreamContour, 'k.', xAeroStreamContour,...
%rAeroStreamContour, 'y.', xAeroStreamContourTrunc,...
%rAeroStreamContourTrunc, 'c.', xAeroExpansion,...
%rAeroExpansion, 'y.', xAeroStream, rAeroStream, 'r.',...
%xAeroStreamCowl, rAeroStreamCowl, 'r.', xAeroStreamTrunc,...
%rAeroStreamTrunc, 'b.');
%else
%plot(xSADcontour, rSADcontour, 'g.', xStreamContour,...
%rStreamContour, 'k.', xAeroStreamContour,...
%rAeroStreamContour, 'y.', xAeroExpansion, rAeroExpansion,...
%'y.', xAeroStream, rAeroStream, 'r.', xAeroStreamCowl,...
%rAeroStreamCowl, 'r.');
%end

%disp('MaSAD(PointNum)');
%disp(MaSAD(PointNum));
%disp(' ');
%disp('rAeroExpansion');
%disp(rAeroExpansion);
%disp(' ');
%disp('Area Ratio');
%disp((rStreamContour(ll))^2);
%disp(' ');