
clear;
%******************************************************************************************************
%                          PROBLEM PARAMETER INPUTS                                                   %                                                                  
%******************************************************************************************************
T_c   = 3216.96;            % Temperature in the combustion chamber (K)
P_c   = 800*6894.757;   % Pressure in the combustion chamber (Pa) 
altitude = 0;         % Target altitude (meter)
h_th  = 0.02413;          % Throat  diameter (meter) == 0.846 inches

% The constants to calculate the specific hear gamma 
a1 =  0.144236559804932;
a2 = -2.92510283447229e-05;
a3 =  0.32204378095657e-06;
a4 = -1.71921724162008e-10;
a5 =  4.43350498618313e-14;
a6 = -4.29053413942424e-18;

% Method of Characteristics
num     = 10;            % Number of Characteristic lines
theta_i = 5;             % Initial step in theta (deg) 
R_ex    = 2;             % Expansion radius (R/Rt). 
teta_N_init = 20;        % Expansion angle (deg)                         

% Mach number at the exit
Select       = 1;        % Set to '1' to calculate the ideal exit mach number
Desired_mach = 3.35;     % If Select is not set to '1', enter the desired exit mach number

%Pintle Option
pintle_op    = 2;        % Set to '1' if there is a pintle, set to '2' for no pintle

%Plot
plotter1 = 1;            % Set to '1' to plot nozzle contour and the Mach lines
plotter2 = 0;            % Set to '1' to plot the thrust curve vs. exit area
plotter3 = 0;            % Set to '1' to plot Mach number and pressure ratio across the nozzle length

%*******************************************************************************************************


%******************************************************************************************************
%                                      PINTLE SHAPES                                                  %    
%******************************************************************************************************

if pintle_op == 2
    pintle_rad = 0.0;    % Pintle tip radius
end    
    
if pintle_op == 1
    fprintf('Choose one of the following pintle geometries')
    fprintf('\n 1. Hemisphere')
    fprintf('\n 2. Elliptical')
    fprintf('\n 3. Conical')
    fprintf('\n 4. Parabolic')
    geometry=input('\nEnter the number of your selection: ');

    %Hemisphere
    if geometry == 1
        RR = input('\nEnter the radius: ');
        LL = RR;
    end                                                                                                                                                                                        

    %Elliptical, Conical, Tangent ogive and Secant ogive
    if geometry > 1 & geometry < 4
        RR = input('\nEnter the height: ');
        LL = input('\nEnter the length: ');
    end

    %Parabolic
    if geometry == 4
        RR = input('\nEnter the height: ');
        LL = input('\nEnter the length: ');
        K = input('\nEnter the value of K (0.5, 0.75 and 1 for 1/2, 3/4 and a full parabola respectively: ');   
    end
    
    pintle_rad = RR;
    
end


%dh_th = [0 0 0 0 0 0 0];
W     = 1;                                                   % Molecular weight (kg/kmol)
dh = h_th/100;
max_iter = 10000;
R = 8314/W;
H_th = h_th;
r_th = h_th/2;
P_amb = (4.35055e-4 * (altitude^2)) - (1.16037e1 * altitude) + 1.01084e5; %N/m2
P_amb = 3000000;
T_amb = 300;                                                 % Ambient temperature (K)


%******************************************************************************************************
%                                   THRUST CALCULATION                                                %                                              
%******************************************************************************************************
%find where P becomes u
teta = T_c / 1000;
Cp_mol = a1 + (a2*teta) + (a3*(teta^2)) + (a4*(teta^3)) + (a5*(teta^4)) + (a6*(teta^5)); % KJ/kmol-K
gamma = Cp_mol / (Cp_mol - 8.315);
width = 0.06;      % Nozzle width (meters)
h(1) = h_th;
A_star = h_th*width;
M =1;
dM1 = .1;
for i=1: max_iter
    h(i) = h(1) + (i-1)*dh;
    Ae(i) = h(i)*width;
    A_Asq = (Ae(i)/A_star)^2;
    A_ratio(i)=sqrt(A_Asq);
    
    %Newton Rhapson on Eq. 5.20 - Anderson text
    res = 1;
    if i > 1
        M = Ma(i-1);
    end
    
     while res > .001
        M2 = M + dM1;
        funa1 = -A_Asq + (1/M^2)*((2/(gamma+1))*(1+(gamma-1)*M^2/2))^((gamma+1)/(gamma-1));
        funa2 = -A_Asq + (1/M2^2)*((2/(gamma+1))*(1+(gamma-1)*M2^2/2))^((gamma+1)/(gamma-1));
        dv_dm = (funa2-funa1)/dM1;
        
        M = M - funa1/dv_dm;
        res = abs(funa1);
        
    end
    Ma(i) = M;
    
    % Find Pressure
    P(i) = P_c*(1+(gamma-1)*Ma(i)^2/2)^(-gamma/(gamma-1));
    
    % Find thrust for each point
    Te(i) = T_c/(1+(gamma-1)*Ma(i)^2/2);
    Tt(i) = T_c/(1+(gamma-1)/2);
    Ve(i) = Ma(i)*sqrt(Te(i)*gamma*R);
    Vt(i) = sqrt(Tt(i)*gamma*R);
    rhot(i) = P(i)/(R*Te(i));
    mdot(i) = rhot(i)*Ve(i)*Ae(i);
    TT(i) = mdot(i)*Ve(i) + (P(i) - P_amb)*Ae(i);
    
    teta = Te(i) / 1000;
    Cp_mol = a1 + (a2*teta) + (a3*(teta^2)) + (a4*(teta^3)) + (a5*(teta^4)) + (a6*(teta^5)); % KJ/kmol-K
    gamma = Cp_mol / (Cp_mol - 8.315);
    
    if P(i) < P_amb
        %break
        %Calculate the pressure if shock wave exists at the exit plane
        P_exit = P(i)*(1+(gamma*2/(gamma+1))*(Ma(i)^2-1));
        
         if P_exit <= P_amb
             P(i) = P_exit;
             break
         else
         end
        
    else
    end
    
end

if plotter2 == 1
    figure(2)
    plot(Ae,TT)
    title('Thrust curve')
    xlabel('Exit Area (m^2)')
    ylabel('Thrust (N)')
end
%******************************************************************************************************


%******************************************************************************************************
%                                 MAX THRUST MIN LENGTH                                               %                                              
%******************************************************************************************************
%   Determine the nominal exit area of the nozzle 
%   to maximize thrust

[a,b]=max(TT);
%   Over or Underexpand the nozzle
b = b;
A_max = Ae(b);
Max_thrust = TT(b);

if plotter2 == 1
hold on;
plot(A_max,Max_thrust,'r*')
legend('Thrust Curve','Max Thrust')
end

%******************************************************************************************************



%******************************************************************************************************
%                                 METHOD OF CHARACTERISTICS                                           %                                              
%******************************************************************************************************

if Select == 1
    M_e = Ma(b);     %Mach number at ideal exit 
    else
    M_e = Desired_mach;
end

%Find theta_max by using equation 11.33
theta_max = (180/pi)*(sqrt((gamma+1)/(gamma-1))*atan((sqrt((gamma-1)*(M_e^2-1)/(gamma+1))))-atan(sqrt(M_e^2-1)))/2;

%  D_theta for each char line
del_theta = (theta_max - theta_i)/(num-1);

% Find 
T = T_c;
for i=1:num
    %   Initialize mach numeber

    for j=1:num
        if i==1
            %Theta for each line (first lines)
            theta(i,j) = theta_i + del_theta*(j-1);
            nu(i,j) = theta(i,j);
            K_m(i,j) = theta(i,j) + nu(i,j);
            K_p(i,j) = theta(i,j) - nu(i,j);       
        
        elseif i > 1
        
            K_p(i,j) = -K_m(1,i);
            
            % Find Thetas
            if j >= i
                theta(i,j) = del_theta*(j-i);
            else 
                %theta(i,j) = theta(j,i-1);
                theta(i,j) = theta(j,i);
                
            end
            nu(i,j) = theta(i,j) - K_p(i,j);
            K_m(i,j) = theta(i,j) + nu(i,j);
        end
    
    % Prandtl-Meyer function (using Newton Rhapson
    dM = .1; % Leave at about .1
    if j == 1
        M_ex(i,j) = 1.00; 
        %   Constant Specific Heat
        T(i,j) = T_c; 
        teta(i,j) = T(i,j) / 1000;
        Cp_mol(i,j) = a1 + (a2*teta(i,j)) + (a3*(teta(i,j)^2)) + (a4*(teta(i,j)^3)) + (a5*(teta(i,j)^4)) + (a6*(teta(i,j)^5)); % KJ/kmol-K
        %Cp(i,j) = a1 + (a2 * T(i,j)) + (a3 * (T(i,j)^2)) + (a4 * (T(i,j)^3)) + (a5 * (T(i,j)^4));  % J/kg-K
        gamma(i,j) = Cp_mol(i,j) / (Cp_mol(i,j) - 8.315);
    else
        M_ex(i,j) = M_ex(i,j-1);
        gamma(i,j) = gamma(i,j-1);
    end
        M(i,j) = M_ex(i,j);    

        
    res = 1;
    while res > .01
        M2(i,j) = M(i,j) + dM;
        funv1(i,j) = (-nu(i,j)*(pi/180)+(sqrt((gamma(i,j)+1)/(gamma(i,j)-1))*atan((sqrt((gamma(i,j)-1)*(M(i,j)^2-1)/(gamma(i,j)+1))))-atan(sqrt(M(i,j)^2-1))));
        funv2(i,j) = (-nu(i,j)*(pi/180)+(sqrt((gamma(i,j)+1)/(gamma(i,j)-1))*atan((sqrt((gamma(i,j)-1)*(M2(i,j)^2-1)/(gamma(i,j)+1))))-atan(sqrt(M2(i,j)^2-1))));
        dv_dm(i,j) = (funv2(i,j)-funv1(i,j))/dM;
        
        M(i,j) = M(i,j) - funv1(i,j)/dv_dm(i,j);
        res(i,j) = abs(funv1(i,j));
        
    end
        M_ex(i,j) = M(i,j);  
    
    if j ~= 1
        %   Constant Specific Heat
        T(i,j) = T(i,j-1) * (1 + (((gamma(i,j-1) - 1)/2) * (M_ex(i,j-1) ^ 2))) / (1 + (((gamma(i,j-1) - 1)/2) * (M_ex(i,j) ^ 2))); 
        teta(i,j) = T(i,j) / 1000;
        Cp_mol(i,j) = 4.184 * (a1 + (a2*teta(i,j)) + (a3*(teta(i,j)^2)) + (a4*(teta(i,j)^3)) + (a5*(teta(i,j)^-2))); % KJ/kmol-K
        %Cp(i,j) = a1 + (a2 * T(i,j)) + (a3 * (T(i,j)^2)) + (a4 * (T(i,j)^3)) + (a5 * (T(i,j)^4));  % J/kg-K
        gamma(i,j) = Cp_mol(i,j) / (Cp_mol(i,j) - 8.315); 
    end
    
    % Find the angle mu
    mu(i,j) = (180/pi)*asin(1/M_ex(i,j));
         
    end
    
    % Add last point to char line
    theta(i,num+1) = theta(i,num);
    nu(i,num+1) = nu(i,num);
    K_m(i,num+1) = K_m(i,num);
    K_p(i,num+1) = K_p(i,num);
    
end

char = zeros(num,num+1,2);
 
for i = 1 : num
 teta_N(i) = i * teta_N_init ./ num;
end
   
R_ex = R_ex * (h_th/2);


%****************************************** NO PINTLE SECTION **********************************************************************    
    for i=1:num    
    for j=1:num+1
        
% Draw points of intersection            
        %   Point 1 of all char lines          
        if j == 1 
            char(i,j,1) = R_ex * sin(teta_N(i) * pi / 180);
            char(i,j,2) = (h_th/2 + R_ex) - (R_ex * cos(teta_N(i) * pi / 180));
        end

        %   Where first line hits the symmetry line
        if i == 1 & j==2            
            char(i,j,1) = 0.006;
            char(i,j,2) = 0;
        end
        
        %   Where all other lines hit the symmetry line
        if j == i+1 & j > 2           
              char(i,j,1) = -char(i-1,j,2)/tan((pi/180)*(.5*theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)))) + char(i-1,j,1);
              char(i,j,2) = 0;
              test(i,j) = (theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)));
              testpty(i,j) = char(i-1,j,2);
              testptx(i,j) = char(i-1,j,1);                
        end
            
        %   All other data points for char 1 calculated
        if i ==1 & j>2 & j ~= i+1
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(j-1,1)+theta(i,j-1))-.5*(mu(j-1,1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(1,1,2) - char(1,1,1)*C_m;
                char(1,j-1,2) - char(1,j-1,1)*C_p];
                iterm(:,1)=inv(A)*B;
                char(i,j,1) = iterm(2,1);
                char(i,j,2) = iterm(1,1);
        end
        
        %   All other points for all char lines calculated
        if i > 1 & j~=i+1 & j>2        
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(i-1,j-1)+theta(i,j-1))-.5*(mu(i-1,j-1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(i-1,j,2) - char(i-1,j,1)*C_m; char(i,j-1,2) - char(i,j-1,1)*C_p];
                
            iterm(:,1) = inv(A)*B;
            char(i,j,1) = iterm(2,1);
            char(i,j,2) = iterm(1,1);  
        end
    end
end
    
%****************************************** WITH PINTLE SECTION **********************************************************************

for i = 1 : num
    h_th(i) = H_th;   % Throat  height (meters) == 0.846 inches
end  

if pintle_op == 1

    r_new = RR;
    for i=1:num
    for j=1:num+1
       if j == i+1 & j>1 
           start=0;
    for x = start:LL/500:LL
        %Hemisphere
        if geometry == 1
            r_new = RR .* sqrt(1 - ((x^2)/(LL^2)));
        end
        %Elliptical
        if geometry == 2
            r_new = RR .* sqrt(1 - ((x^2)/(LL^2)));
        end
        %Conical
        if geometry == 3
            r_new = RR - (x .* RR ./ LL);
        end
        %Tangent ogive
        if geometry == 100
            rho = ((RR^2) + (LL^2)) / (2*RR);
            r_new = RR -(sqrt((rho^2)-((x-LL)^2)) + RR - rho);
        end
        %Secant ogive
        if geometry == 100
            rho = ((RR^2) + (LL^2)) / (2*RR);
            r_new = sqrt((rho^2)-((x-LL)^2)) + RR - rho;
        end
        %Parabolic
        if geometry == 4
            r_new = RR * ( ((2*(x/LL)) - (K*((x/LL)^2))) / (2-K));
        end
        angle = atan((char(i,j,1)/r_th));
        y_new = x / tan(angle);
        tolerance = r_th - (y_new + r_new);
        if tolerance < r_th/100000
            fprintf('Tolerance met');
            break
        end
    XX = x;
    h_th(i) = (r_th - r_new) * 2;
    dh_th1(i+1) = r_new;
    dh_th1(1) = RR;
end
end
end
end



for AA = 1 : num
    for i = AA
        for j = AA : num+2
            dh_th(i,j) = dh_th1(i);
        end
    end
end

for AA = 1 : num
    for j = AA
        for i = AA : num
            dh_th(i,j) = dh_th1(j);
        end
    end
end
 
dh_th(num,num) = dh_th(num,num-1);
dh_th(num-1,num+2) = 0;

for i=1:num           
    for j=1:num+1
        
% Draw points of intersection    

        %   Point 1 of all char lines          
        if j == 1 
            char(i,j,1) = 0;
            char(i,j,2) = r_th - RR;
        end
        
        %   Where first line hits the symmetry line
        if i == 1 & j==2           

            char(i,j,1) = (-h_th(i)/2)/tan((pi/180)*(theta(1,j-1)-mu(1,j-1)));
            char(i,j,2) = 0; 
        end
       
        %   Where all other lines hit the symmetry line
        if j == i+1 & j>2       
            
              char(i,j,1) = -char(i-1,j,2)/tan((pi/180)*(.5*theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)))) + char(i-1,j,1);
              char(i,j,2) = 0;
              test(i,j) = (theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)));
              testpty(i,j) = char(i-1,j,2);
              testptx(i,j) = char(i-1,j,1);                
        end
        
        %   All other data points for char 1 calculated
        if i ==1 & j>2 & j ~= i+1
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(j-1,1)+theta(i,j-1))-.5*(mu(j-1,1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(1,1,2) - char(1,1,1)*C_m;
                char(1,j-1,2) - char(1,j-1,1)*C_p];
                iterm(:,1)=inv(A)*B;
                char(i,j,1) = iterm(2,1);
                char(i,j,2) = iterm(1,1);
        end
        
        %   All other points for all char lines calculated
        if i > 1 & j~=i+1 & j>2        
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(i-1,j-1)+theta(i,j-1))-.5*(mu(i-1,j-1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(i-1,j,2) - char(i-1,j,1)*C_m; char(i,j-1,2) - char(i,j-1,1)*C_p];
                
            iterm(:,1) = inv(A)*B;
            char(i,j,1) = iterm(2,1);
            char(i,j,2) = iterm(1,1);  
        end        
    end
end
end

%  Fill in similar points (where char lines share points)
for i = 2:num
    for j=2:num
        char(j,i,1) = char(i-1,j+1,1);
        char(j,i,2) = char(i-1,j+1,2);
    end
end
        
% *******************************Make the nozzle shape and extend the char lines to wall********************

%   Initial start point of the nozzle (at throat)
noz(1,1) = R_ex * sin(teta_N(i) * pi / 180);
noz(1,2) = (h_th(i)/2 + R_ex) - (R_ex * cos(teta_N(i) * pi / 180));

%   Find all the points of the nozzle
for i = 2 : num
    %   Find different slopes and points to intersect
    m1 = tan((pi/180)*(theta(i-1,num)+mu(i-1,num)));    
    if i ==2
        m2 = (pi/180)*theta_max;
    else
        m2 = ((pi/180)*(theta(i-1,num+1)));
    end
    m3 = ((pi/180)*(theta(i-1,num)));
    m4 = tan((m2+m3)/2);
    
    A = [1,-m4; 1,-m1];
    B = [noz(i-1,2) - noz(i-1,1)*m4; char(i-1,num+1,2) - char(i-1,num+1,1)*m1];
                
    iterm(:,1) = inv(A)*B;
    noz(i,1) = iterm(2,1);
    noz(i,2) = iterm(1,1); 
    
    %   Extend char lines to wall
    char(i-1,num+2,1)= noz(i,1);
    char(i-1,num+2,2)= noz(i,2);
end

%Last line
m1 = tan((pi/180)*(theta(num,num)+ mu(num,num)));
m2 = ((pi/180)*(theta(num-1,num)));
m3 = ((pi/180)*(theta(num,num+1)));

m4 = tan((m2+m3)/2);
A = [1,-m4; 1,-m1];
B = [noz(num,2) - noz(num,1)*m4; char(num,num+1,2) - char(num,num+1,1)*m1];
                
iterm(:,1) = inv(A)*B;
noz(num+1,1) = iterm(2,1);
noz(num+1,2) = iterm(1,1); 

    
%   Extend char lines to wall
char(num,num+2,1)= noz(num+1,1);
char(num,num+2,2)= noz(num+1,2);

%******************************************************************************************************
%                                       PLOTTING SECTION                                              %                                              
%******************************************************************************************************  
if pintle_op == 1   
    ddh_th = dh_th1;
    ddh_th(1) = 0;
    ddh_th(num-1) = 0;
end

if plotter1 == 1
%   Plot for loop for char lines
figure(1);
%subplot(2,1,1)
for i = 1 : num
    
    figure(1)
    hold on;
    if pintle_op == 2
        plot(char(i,:,1)/r_th,char(i,:,2)/r_th,'-')
    end
    if pintle_op == 1
        plot(char(i,:,1),dh_th(i,:)+char(i,:,2),'-')
    end
    hold on;
    if pintle_op == 2
        plot(char(i,:,1)/r_th,-char(i,:,2)/r_th,'-')
    end
    if pintle_op == 1   
        plot(char(i,:,1),-dh_th(i,:)-char(i,:,2),'-')
    end    
end

if pintle_op == 1
figure(1)
for x = 0:LL/1000:LL
    %Hemisphere
    if geometry == 1   
        y = RR .* sqrt(1 - ((x^2)/(LL^2)));
    end
    %Elliptical
    if geometry == 2
        y = RR .* sqrt(1 - ((x^2)/(LL^2)));
    end
    %Conical
    if geometry == 3
        y = RR - (x .* RR ./ LL);
    end
    %Tangent ogive
    if geometry == 100
        rho = ((RR^2) + (LL^2)) / (2*RR);
        y = RR -(sqrt((rho^2)-((x-LL)^2)) + RR - rho);
    end
    %Secant ogive
    if geometry == 100
        rho = ((RR^2) + (LL^2)) / (2*RR);
        y = sqrt((rho^2)-((x-LL)^2)) + RR - rho;
    end
    %Parabolic
    if geometry == 4
        y = RR * ( ((2*(x/LL)) - (K*((x/LL)^2))) / (2-K));
    end
    plot(x,y,'r-')
    hold on
    plot(x,-y,'r-')
    hold on
end
end   

%   Plot the nozzle shape
figure(1)
%subplot(1,2,1)
hold on;
if pintle_op == 2
    plot(noz(:,1)/r_th,noz(:,2)/r_th,'k')
end
if pintle_op == 1
    plot(noz(:,1),ddh_th(:)+noz(:,2),'k')
end
hold on;
if pintle_op == 2
plot(noz(:,1)/r_th,-noz(:,2)/r_th,'k')
end
if pintle_op == 1
    plot(noz(:,1),-ddh_th(:)-noz(:,2),'k')
end
hold on;
teta_N = 0;
for i = 1 : teta_N_init/0.01    
    nozx(i) = R_ex * sin(teta_N * pi / 180);
    nozy(i) = (H_th/2 + R_ex) - (R_ex * cos(teta_N * pi / 180));
    teta_N = teta_N + 0.01;
end
plot(nozx/r_th,nozy/r_th,'k')
hold on
plot(nozx/r_th,-nozy/r_th,'k')
hold on
    
[a,b] = max(noz);
%plot(a(1),A_max/width/2,'g*')
hold on;
%plot(a(1),-A_max/width/2,'g*')
%title('Max Thrust (minimum length) Nozzle Design')
title('Nozzle Contour - R/Rt = 2.0')
xlabel('Nozzle length ')
ylabel('Nozzle height ')
grid on
%legend('Nozzle shape','Area_e_x_i_t(predicted)')
 
end


%   Find  % errors in A/A* and Mexit
error_Area = 100*(width*2*noz(num,2) - A_max)/(A_max)
error_Mach = 100*(M_e - M_ex(num,num))/M_e

teta = T_c / 1000;
Cp_mol = a1 + (a2*teta) + (a3*(teta^2)) + (a4*(teta^3)) + (a5*(teta^4)) + (a6*(teta^5)); % KJ/kmol-K
gamma = Cp_mol / (Cp_mol - 8.315);


%   Plot Mach Number and pressure through nozzle
Mnoz(1) = 1.0;  %   Choked Flow
M = Mnoz(1);
for i=1: size(noz,1)
    Ae(i) = 2*noz(i,2)*width;
    A_Asq = (Ae(i)/A_star)^2;
    A_ratio(i)=sqrt(A_Asq);
    
    %Newton Rhapson on Eq. 5.20 - Anderson text
    res = 1;
    if i > 1
        M = Mnoz(i-1);
    iter=0;
    
     while res > .001
        iter = iter + 1; 
        M2 = M + dM1;
        funa1 = -A_Asq + (1/M^2)*((2/(gamma+1))*(1+(gamma-1)*M^2/2))^((gamma+1)/(gamma-1));
        funa2 = -A_Asq + (1/M2^2)*((2/(gamma+1))*(1+(gamma-1)*M2^2/2))^((gamma+1)/(gamma-1));
        dv_dm = (funa2-funa1)/dM1;
        
        M = M - funa1/dv_dm;
        res = abs(funa1);
        if iter > 300
            break
        end
        
        
    end
    Mnoz(i) = M;
    Te   = T_c/(1+(gamma-1)*Mnoz(i)^2/2);
    teta = Te / 1000;
    Cp_mol = a1 + (a2*teta) + (a3*(teta^2)) + (a4*(teta^3)) + (a5*(teta^4)) + (a6*(teta^5)); % KJ/kmol-K
    gamma = Cp_mol / (Cp_mol - 8.315);
    
    end
    % Find Pressure
    Pnoz(i) = P_c*(1+(gamma-1)*Mnoz(i)^2/2)^(-gamma/(gamma-1));
    Pback = Pnoz(i);
end

Pback = Pback

if plotter3 == 1
    figure(3);
    %subplot(1,2,2)
    plot(noz(:,1),Mnoz,'b')
    hold on;
    plot(noz(:,1),Pnoz/P_amb,'r')
    hold on;
    plot(noz(size(noz,1),1),M_e,'k*')
    hold on;
    plot(noz(size(noz,1),1),1,'g*')
    axis([0 0.35 0 5])
    grid on
    legend('Mach Number','P/P_a_m_b','M_e_x_i_t(predicted)','P_a_m_b/P_a_m_b')
end