clc
clear all
format long

%Global Variables%
%starting and ending years (compared to Jan, 2007)
%first year is the real year, subtract 2007 for the math to work
global STARTYEAR
global ENDYEAR
STARTYEAR = 2030-2007;
ENDYEAR = 2040-2007;

%Loop counters, calculated so that each time slice is 1/2 of an Earth Day
global NUM
global BANDS
NUM = (ENDYEAR-STARTYEAR)*365*2;
BANDS = 4;

%CONSTANTS
NODES = 6;
HOPS = 2;



[x,y,t] = posnCalc();
%Two hop space routing. Assumes all Tx's and Rx's are the same.
%Earth is one destination, Mars is the other.
source = 1;
target = 2;

%Calculate the distance in m as a comparison
%Calculate the angles as a comparison
for i = 1:NUM
    dist(i) = sqrt((x(1,i)-x(2,i)).^2 + (y(1,i)-y(2,i)).^2);
    thetaEM(:,i) = anglesCalc([x(1,i),y(1,i)],[x(2,i),y(2,i)]);
    thetaL4M(:,i) = anglesCalc([x(3,i),y(3,i)],[x(2,i),y(2,i)]);
    thetaL5M(:,i) = anglesCalc([x(4,i),y(4,i)],[x(2,i),y(2,i)]);
end

%Use this variable to calculate the nubmer of days on Mars with 0 service.
zeroDays = 0;

distTable = zeros(NODES);
jumps = zeros(3,NODES^HOPS);
oneHop = zeros(2,NODES);
index = zeros(1,NUM);
%Run this once for each time slice
for j = 1:NUM
    %Build the distances table
    for i = 1:NODES
        for k = 1:NODES
            distTable(i,k) = sqrt((x(i,j)-x(k,j)).^2 + (y(i,j)-y(k,j)).^2);
            %Add the effects of the Sn on the direct route
            sunTheta = anglesCalc([x(i,j),y(i,j)],[x(k,j),y(k,j)]);
            sunMult = max(sunEffect(sunTheta(1)),sunEffect(sunTheta(2)));
            distTable(i,k) = distTable(i,k)*sunMult;
        end
    end
    %Find the direct route
    direct = distTable(source,target);
    Rb(1,j) = linkBudget(direct);
    if(Rb(1,j)<=1)
        zeroDays = zeroDays + 1;
    end
    
    
    %Find the best route with 1 hop and only Earth L4
    for i = 1:3
        %Find all the first hops from the source
        earthL4Hop(1,i) = distTable(source,i); 
        %Find all the distances from the first hop to the target
        earthL4Hop(2,i) = distTable(i, target);
    end
    earthL4HopMin = min(max(earthL4Hop));
    Rb(2,j) = linkBudget(earthL4HopMin);
    
    %Find the best route with 1 hop and only Earth L5
    for i = 1:3
        %Introduce a secondary variable so the loop can pick L5, not L4
        %1 = Earth, 2 = Mars, 4 = L5
        z = [1,2,4];
        %Find all the first hops from the source
        earthL5Hop(1,i) = distTable(source,z(i)); 
        %Find all the distances from the first hop to the target
        earthL5Hop(2,i) = distTable(z(i), target);
    end
    earthL5HopMin = min(max(earthL5Hop));
    Rb(3,j) = linkBudget(earthL5HopMin);
    
      
    %Find the best route with 1 hop and only Earth Lagranges
    for i = 1:4
        %Find all the first hops from the source
        earthHop(1,i) = distTable(source,i); 
        %Find all the distances from the first hop to the target
        earthHop(2,i) = distTable(i, target);
    end
    earthHopMin = min(max(earthHop));
    Rb(4,j) = linkBudget(earthHopMin);
    
    %Find the best route with 1 hop
    for i = 1:NODES
        %Find all the first hops from the source
        oneHop(1,i) = distTable(source,i); 
        %Find all the distances from the first hop to the target
        oneHop(2,i) = distTable(i, target);
    end
    oneHopMin = min(max(oneHop));
    Rb(5,j) = linkBudget(oneHopMin);
    
    %Find the best path with 2 hops
    for i = 1:NODES
        for k = 1:NODES
        %Find all the first hops from the source
        jumps(1,(i-1)*NODES+k) = distTable(source,i); 
        %Find all the  second hops from that first hop
        jumps(2,(i-1)*NODES+k) = distTable(i,k);
        %Find all the distances from the second hop to the target
        jumps(3,(i-1)*NODES+k) = distTable(k, target);
        end
    end
    %Find the longest single leg of each path
    jumpMax = max(jumps);
    %Find the path with the shortest maximum leg distance
    [jumpMin, index(j)] = min(jumpMax);
    %Find the max bandwidth of that path
    Rb(6,j) = linkBudget(jumpMin);
end
%Plug it into Link budget to find bandwidth.
%Rb = linkBudget(dist);
%Plot Things
%Convert time from seconds to Earth Years%
t = t./(365*24*3600)+2007;
%Average Data Rate (Mbps)
for i = 1:5
    avgRb(i) = sum(Rb(i,:))/NUM;
    minRb(i) = min(Rb(i,:));
end

%Calculate the increase to average and minimum bandwidth, in percentage.
for i = 1:4
    incAvg(i) = ((avgRb(1+i)/avgRb(1))*100)-100;
    incMin(i) = ((minRb(1+i)/minRb(1))*100)-100;
end

figure('Name','Data Link Rates');
subplot(5,1,1)
plot(t,dist/10^3,'r');
title('Earth-Mars Distance (In Km)');
subplot(5,1,2)
plot(t,Rb(1,:),'c');
subplot(5,1,3)
plot(t,Rb(2,:),'g');
title('One Hop with just Earth L4 available (In Mbps)');
subplot(5,1,4)
plot(t,Rb(4,:),'b');
title('One Hop with just Earth L5 available  (In Mbps)');

title('Direct Link to Earth (In Mbps)');
subplot(5,1,5)
plot(t,Rb(4,:),'m');
title('One Hop with Earth L4 and L5 Satellites (In Mbps)');

% subplot(5,1,4)
% plot(t,oneHopRb(4,:),'g');
% title('One Hop with Earth and Mars L4 and L5 Satellites (In Mbps)');
% subplot(5,1,5)
% plot(t,twoHopRb(4,:),'b');
% title('Two Hops with Earth and Mars L4 and L5 Satellites (In Mbps)');

figure('Name','Sun Effects');
subplot(4,1,1)
plot(t,thetaEM(1,:),'b');
title('Look Angle between the Sun and Mars from Earth (deg)');
subplot(4,1,2)
plot(t,thetaEM(2,:),'r');
title('Look Angle between the Sun and Earth from Mars (deg)');
subplot(4,1,3)
plot(t,thetaL4M(1,:),'g');
title('Look Angle between the Sun and Mars from Earth L4 (deg)');
subplot(4,1,4)
plot(t,thetaL5M(1,:),'r');
title('Look Angle between the Sun and Mars from Earth L5 (deg)');



function [x,y,t] = posnCalc()

%CONSTANTS%
%number of iterations (unitless)
global NUM
%time range
global STARTYEAR
global ENDYEAR
%GIVENS%
%Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/%
%Start time is 2007 Jan 15, 0030 UT, Earth Time of Periapsis
%Earth is 1, Mars is 2
%Perihelion (m)
r_p = [147.10*10^9, 206.62*10^9];
%Aphelion (m)
r_a = [152.10*10^9, 249.23*10^9];
%Standard Gravitational Parameter of Sun (m^3/s^2)
mu = 1.327124400*10^20;
%Time Range (seconds)
%Counting up from 2007 Jan 15 at 0030 UT.
t = linspace(STARTYEAR*365*24*3600,ENDYEAR*365*24*3600,NUM);

%ORBIT CHARACTERISTICS - from Curtis, Table 3.1%
%Eccentricity (unitless)
ecc = (r_a-r_p)./(r_a+r_p);
%semimajor & semiminor axis (m)
a = (r_p+r_a)./2;
b = a.*sqrt(1-ecc.^2);
%Angular Momentum (m^2/s)
h = sqrt(r_p.*(1+ecc)*mu);
%period (seconds)
T = (((2*pi)/(mu^2)).*(h./(sqrt(1-ecc.^2))).^3);

%VARIABLES%
%True Anomaly (radians)
theta_base = linspace(0,2*pi,NUM);

%CALCULATIONS%
%180deg offset so it starts at perihelion when t = 0.
t_new(1,:) = T(1)/2+t;
%Martian perihelion occurs on 01 Jun 2007 at 0720
%This is 137 days, 6 hours, and 50 minutes after Earth perihelion
t_mars = 137*24*3600 + 6*3600 + 50*60;
t_new(2,:) = T(2)/2+t-t_mars;
%Mean Anomaly-eq3.12 (rad)
for i = 1:2
    M_e(i,:) = 2*pi*t_new(i,:)/T(i);
end
%Solve Keplers Equation%
%Initial guess for E
E = zeros(2,NUM);
for i = 1:2
    if(M_e(i) > pi)
        E(i,:) = M_e(i,:)-ecc(i)/2;
    else
        E(i,:) = M_e(i,:)+ecc(i)/2;
    end
end

% Set the acceptable error
error = 10^(-10);
diff = [1,1];

%iterate until acceptable error is reached
for i = 1:2
    while (diff(i)>error)
        diff(i) = ((E(i)-ecc(i)*sin(E(i))-M_e(i))./(1-ecc(i)*cos(E(i))));
        E(i) = E(i) - diff(i);
    end
end

%Argument of Periapsis + RAAN (rad)
alpha = [102.9*pi/180,(286.5+49.6)*pi/180];
%find the true anomaly
theta = 2*atan((sqrt((1+ecc)/(1-ecc))*tan(E/2)));
%Offset from Solar Center
for i = 1:2
    C_x(i) = (a(i) - r_p(i))*cos(alpha(i));
    C_y(i) = (a(i) - r_p(i))*sin(alpha(i));
end

%Determine complete orbital ellipse
for j = 1:2
    for i = 1:NUM
        x(j,i) = a(j)*cos(alpha(j))*cos(theta_base(i))-b(j)*sin(alpha(j))*sin(theta_base(i))+C_x(j);
        y(j,i) = a(j)*sin(alpha(j))*cos(theta_base(i))+b(j)*cos(alpha(j))*sin(theta_base(i))+C_y(j);
    end
end

%Earth's position at given time
x(1,:) = a(1)*cos(alpha(1))*cos(theta(1,:))-b(1)*sin(alpha(1))*sin(theta(1,:))+C_x(1);
y(1,:) = a(1)*sin(alpha(1))*cos(theta(1,:))+b(1)*cos(alpha(1))*sin(theta(1,:))+C_y(1);
%Mars' position at given time
x(2,:) = a(2)*cos(alpha(2))*cos(theta(2,:))-b(2)*sin(alpha(2))*sin(theta(2,:))+C_x(2);
y(2,:) = a(2)*sin(alpha(2))*cos(theta(2,:))+b(2)*cos(alpha(2))*sin(theta(2,:))+C_y(2);
%Roughly Earth's L4 position at given time%
%TRUE ANOMALY OFFSET INSTEAD OF MEAN - ERROR INTRODUCED%
L4 = 60*pi/180;
x(3,:) = a(1)*cos(alpha(1))*cos(theta(1)+L4)-b(1)*sin(alpha(1))*sin(theta(1)+L4)+C_x(1);
y(3,:) = a(1)*sin(alpha(1))*cos(theta(1)+L4)+b(1)*cos(alpha(1))*sin(theta(1)+L4)+C_y(1);
%Roughly Earth's L5 position at given time%
%TRUE ANOMALY OFFSET INSTEAD OF MEAN - ERROR INTRODUCED%
L5 = 60*pi/180;
x(4,:) = a(1)*cos(alpha(1))*cos(theta(1)-L5)-b(1)*sin(alpha(1))*sin(theta(1)-L5)+C_x(1);
y(4,:) = a(1)*sin(alpha(1))*cos(theta(1)-L5)+b(1)*cos(alpha(1))*sin(theta(1)-L5)+C_y(1);
%Roughly Mars' L4 position at given time%
%TRUE ANOMALY OFFSET INSTEAD OF MEAN - ERROR INTRODUCED%
x(5,:) = a(2)*cos(alpha(2))*cos(theta(2)+L4)-b(2)*sin(alpha(2))*sin(theta(2)+L4)+C_x(2);
y(5,:) = a(2)*sin(alpha(2))*cos(theta(2)+L4)+b(2)*cos(alpha(2))*sin(theta(2)+L4)+C_y(2);
%Roughly Mars' L5 position at given time%
%TRUE ANOMALY OFFSET INSTEAD OF MEAN - ERROR INTRODUCED%
x(6,:) = a(2)*cos(alpha(2))*cos(theta(2)-L5)-b(2)*sin(alpha(2))*sin(theta(2)-L5)+C_x(2);
y(6,:) = a(2)*sin(alpha(2))*cos(theta(2)-L5)+b(2)*cos(alpha(2))*sin(theta(2)-L5)+C_y(2);
end

function Rb = linkBudget(dist)
%INFO%
%Dist is an array of NUM values for distance between 2 points%
%Calculates max theoretical data rate at some distance d%
%1 is X, 2 is Ka, 3 is EHF-low, 4 is EHF-high%

%CONSTANTS%
%Number of bands used (unitless)
global BANDS
%Number of iterations on distance
global NUM
%Speed of Light (m/s)
c = 3*10^8;
%Boltsmann's constant (dB(W/k/Hz))
k = -228.5991672;

%GIVENS%
%Frequency (MHz)
f = [8400,32000,1*10^5, 2*10^5];
%Wavelength (m)
lambda = c./(f.*10^6);
%Power (Watts)
P = [5000, 5000,5000,5000];
%Power (dBW)
P_dB = 10*log10(P);
%Transmitter Line Loss (dB, unitless)
LL = 1;
%Transmitter Antenna Diameter (m)
d_Tx = [8.5,8.5,8.5,8.5];
%Transmitter and Receiver Antenna Efficiency (unitless)
eff = 0.65;
%Transmitter Antenna Gain (dB)
Gt = 10*log10(((pi.*d_Tx./lambda).^2).*eff);
%Pointing Loss (dB, unitless)
PL = 1;
%Free Space Loss (dB, unitless)

for j = 1:BANDS
    %note, distance needs to be converted to Km%
    FSL(j) = 32.4+20*log10(dist/(1*10^3))+20*log10(f(j));
end

%Receiver Antenna Diameter (m)
d_Rx = [8.5,8.5,8.5,8.5];
%Receiver Antenna Gain (dB)
Gr = 10*log10(((pi.*d_Rx./lambda).^2).*eff);
%System Noise Temperature (K) - ASSUMPTION
T = [50,80,100,120];
%Sytem Noise Temp (dBK)
T_dB = k+10*log10(T);
%Required Eb/N0 (dB0
EbN0_req = 1;
% Eb/N0 Margin (dB)
EbN0_margin = 3;

%Calculate Data Rate
PtN0 = zeros(1,BANDS);
EbN0 = zeros(1,BANDS);

%Update this to just consider Ka
%for j = 1:BANDS
j = 4;
    %Recieved Pt/N0 (dB-Hz)
    PtN0(j) = P_dB(j) - LL + Gt(j) - PL - FSL(j) + Gr(j) - T_dB(j);
    %Available Eb/N0 (dB)
    EbN0(j) = PtN0(j) - EbN0_req - EbN0_margin;
    %Maximum Theoretical Data Rate (Mb/s)
    Rb = (10^(EbN0(j)/10))/10^6;


end

%This function calculates the three internal angles in a triangle made by:
%The Sun at [0,0], Celestial Body 1 at p1, and Celestial Body 2 at p2
%This is used to determine the nous caused by the Sun
function [theta] = anglesCalc(p1,p2)
    %First calculate the angle between p1 and p2 as seen from the Sun
    theta(3) = acos((dot(p1,p2))/(sqrt(p1(1)^2+p1(2)^2)*sqrt(p2(1)^2+p2(2)^2)))*180/pi;

    %Second calculate the angle between p2 and the Sun, as seen from p1
    % First, center the vector from p1 to p2 on zero
    p1p2 = p2-p1;
    % Second, center the vector from p1 to the Sun on zero
    p1Sun = [0,0]-p1;
    % Then calculate the angle between p1Sun and p1p2
    theta(1) = acos((dot(p1p2,p1Sun))/(sqrt(p1p2(1)^2+p1p2(2)^2)*sqrt(p1Sun(1)^2+p1Sun(2)^2)))*180/pi;
    
    %Finally calculate the angle between p1 and the Sun, as seen from p2
    % First, center the vector from p2 to p1 on zero
    p2p1 = p1-p2;
    % Second, center the vector from p2 to the Sun on zero
    p2Sun = [0,0]-p2;
    % Then calculate the angle between p1Sun and p1p2
    theta(2) = acos((dot(p2p1,p2Sun))/(sqrt(p2p1(1)^2+p2p1(2)^2)*sqrt(p2Sun(1)^2+p2Sun(2)^2)))*180/pi;
    
    %Check to see if all three angles add up to 180
    sum(theta);

end

%This function calculates the additional System Noise Temp (in K) from the Sun
%at any given angle, in degrees. From Eq 8.2 in the DSN
%Telecommunications Link Design Handbook (Slobin, 2015)
%Note that only Ka band is being considered
function [sunMult] = sunEffect(theta)
    % System Noise Temperature (K)
    if (theta <= 0.35)
        T = 10^15;
    elseif (theta <= 0.75)
        T = 1400*exp(-5.1*theta);
    elseif(theta <= 4)
        T = 86*exp(-1.4*theta);
    else
        T = 0;
    end 
    
    %Default System Noise Temp in Ka Band (K)
    T0 = 120;
    %Determine the effect of the Sun's addition to the temp in DB
    sun_dB = 10*log10((T+T0)/T0);
    %Determine the miltiplier for Distance to account for the Sun
    sunMult = 10^(sun_dB/10);
    %Note, going through dB's is uneccessary
    %sunMult = (T+T0)/T0
end



