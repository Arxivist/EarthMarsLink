clc
clear all
format long

%Global Variables%
global NUM
global BANDS
NUM = 10000;
BANDS = 4;

%CONSTANTS
NODES = 6;
HOPS = 2;



[x,y,t] = posnCalc();
%Two hop space routing. Assumes all Tx's and Rx's are the same.
%Earth is one destination, Mars is the other.
source = 1;
target = 2;

%Calculate the distance as a comparison
for i = 1:NUM
    dist(i) = sqrt((x(1,i)-x(2,i)).^2 + (y(1,i)-y(2,i)).^2);
end

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
        end
    end
    %Find the direct route
    direct = distTable(source,target);
    directRb(:,j) = linkBudget(direct);
    
    %Find the best route with 1 hop and only Earth Lagranges
    for i = 1:4
        %Find all the first hops from the source
        earthHop(1,i) = distTable(source,i); 
        %Find all the distances from the first hop to the target
        earthHop(2,i) = distTable(i, target);
    end
    earthHopMin = min(max(earthHop));
    earthHopRb(:,j) = linkBudget(earthHopMin);
    
    %Find the best route with 1 hop
    for i = 1:NODES
        %Find all the first hops from the source
        oneHop(1,i) = distTable(source,i); 
        %Find all the distances from the first hop to the target
        oneHop(2,i) = distTable(i, target);
    end
    oneHopMin = min(max(oneHop));
    oneHopRb(:,j) = linkBudget(oneHopMin);
    
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
    twoHopRb(:,j) = linkBudget(jumpMin);
end
%Plug it into Link budget to find bandwidth.
%Rb = linkBudget(dist);
%Plot Things
%Convert time from seconds to Earth Years%
t = t./(365*24*3600);
%Average Data Rate (Mbps)
avgRb(1) = sum(directRb(4,:))/NUM;
minRb(1) = min(directRb(4,:));
avgRb(2) = sum(earthHopRb(4,:))/NUM;
minRb(2) = min(earthHopRb(4,:));
avgRb(3) = sum(oneHopRb(4,:))/NUM;
minRb(3) = min(oneHopRb(4,:));
avgRb(4) = sum(twoHopRb(4,:))/NUM;
minRb(4) = min(twoHopRb(4,:));

for i = 1:3
    incAvg(i) = ((avgRb(1+i)/avgRb(1))*100)-100;
    incMin(i) = ((minRb(1+i)/minRb(1))*100)-100;
end

figure;
subplot(5,1,1)
plot(t,dist/10^3,'r');
title('Earth-Mars Distance (In Km)');
subplot(5,1,2)
plot(t,directRb(4,:),'c');
title('Direct Link to Earth (In Mbps)');
subplot(5,1,3)
plot(t,earthHopRb(4,:),'m');
title('One Hop with Earth L4 and L5 Satellites (In Mbps)');
subplot(5,1,4)
plot(t,oneHopRb(4,:),'g');
title('One Hop with Earth and Mars L4 and L5 Satellites (In Mbps)');
subplot(5,1,5)
plot(t,twoHopRb(4,:),'b');
title('Two Hops with Earth and Mars L4 and L5 Satellites (In Mbps)');

function [x,y, t] = posnCalc()

%CONSTANTS%
%number of iterations (unitless)
global NUM

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
t = linspace(13*365*24*3600,18*365*24*3600,NUM);

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

%Plot Results
%f = figure;
%plot(t,dist,'b');

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
%Boltsmann's constant (dB(W/k/Hz)
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
Rb = zeros(1,BANDS);
PtN0 = zeros(1,BANDS);
EbN0 = zeros(1,BANDS);
for j = 1:BANDS
    %Recieved Pt/N0 (dB-Hz)
    PtN0(j) = P_dB(j) - LL + Gt(j) - PL - FSL(j) + Gr(j) - T_dB(j);
    %Available Eb/N0 (dB)
    EbN0(j) = PtN0(j) - EbN0_req - EbN0_margin;
    %Maximum Theoretical Data Rate (Mb/s)
    Rb(j) = (10^(EbN0(j)/10))/10^6;
end

end

