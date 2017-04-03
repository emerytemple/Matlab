%**************************************************************************
%* Program will calculate the Mach Number associated with a given
%* Prandtl-Meyer Expansion Angle and Ratio of Specific Heats, Gamma
%* Brandon Denton
%* RIT Graduate Student
%* February 11, 2007
%**************************************************************************
%**************************************************************************
%* Variable Dictionary
%**************************************************************************
%* 1. Mach -------------------- Points Calculated Mach Number
%* 2. MachG ------------------- Guess Mach Number value
%* 3. MaEnd ------------------- Variable that holds the high end value of
%* the calculation range
%* 4. MaStart ----------------- Variable that holds the low end value of
%* the calculation range
%* 5. Mexit ------------------- Desired exit Mach Number of the nozzle
%* 6. vCheck ------------------ Variable that holds the Prandtl-Meyer
%* angle associated with the guess Mach
%* Number
%* 7. vError ------------------ Error and Loop Control Variable
%* 8. vRad -------------------- Acutal Prandtl-Meyer Expansion angle of
%* the current point
%**************************************************************************
%**************************************************************************
%* Start Program
%**************************************************************************
%**************************************************************************
MaStart = MachG;
MaEnd = 100 * Mexit;
vError = 1.0;
while abs(vError) > 1e-10
MachG = (MaEnd + MaStart) / 2;
vCheck = (sqrt((Gamma + 1)/(Gamma - 1)))*...
atan(sqrt(((Gamma - 1)/(Gamma + 1)) *...
((MachG^2) - 1))) - atan((sqrt((MachG^2) - 1)));
vError = vRad - vCheck;
if abs(vError) > 1e-10
if vError > 0.0
MaStart = MachG;
else
MaEnd = MachG;
end
else
Mach = MachG;
end
end