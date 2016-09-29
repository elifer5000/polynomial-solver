function project()
close all; clc; clear;
% Fixed data
A1=[0; 0];
A2=[1; 0];
A3=[0; 1];
side=0.25;

% Input data (0.25..1)
% L1=0.8;
% L2=0.8;
% L3=0.8;
L1=0.6;
L2=0.6;
L3=0.6;
% L1=0.815;
% L2=0.847;
% L3=0.48;
% L1=1;
% L2=0.80;
% L3=0.86;
% L1=0.4330;
% L2=1.0;
% L3=0.5;
% L1=0.8; % 179.999
% L2=0.8;
% L3=0.87423;
% L1=0.6971;
% L2=1.0;
% L3=0.6071;

a= 26278.1+11245.6*L1^4+6125.62*L2^4-17142.6*L3^2+6400.*L3^4+L1^2*(-3514.26-10971.2*L2^2-11520.*L3^2)+L2^2*(-17228.3-1280.*L3^2);
b= -12289.-1499.24*L1^4-2048.*L2^4+5542.56*L3^2+L1^2*(1975.14+3547.24*L2^2-548.76*L3^2)+L2^2*(4939.32+548.76*L3^2);
c= 34895.7+28141.6*L1^4+14829.6*L2^4-32571.3*L3^2+15104.*L3^4+L1^2*(-3372.31-27867.2*L2^2-28416.*L3^2)+L2^2*(-32245.5-1792.*L3^2);
d= -10313.-2998.48*L1^4-4096.*L2^4+7537.89*L3^2+L1^2*(2451.04+7094.48*L2^2-1097.52*L3^2)+L2^2*(7281.89+1097.52*L3^2);
e= 14714.4+22546.4*L1^4+11282.4*L2^4-18761.4*L3^2+11008.*L3^4+L1^2*(1201.4-22820.8*L2^2-22272.*L3^2)+L2^2*(-18401.3+256.*L3^2);
f= -2120.03-1499.24*L1^4-2048.*L2^4+1995.32*L3^2+L1^2*(475.9+3547.24*L2^2-548.76*L3^2)+L2^2*(2342.56+548.76*L3^2);
g= 2000.72+5650.38*L1^4+2578.38*L2^4-3332.68*L3^2+2304.*L3^4+L1^2*(1059.45-5924.76*L2^2-5376.*L3^2)+L2^2*(-3384.12+768.*L3^2);

%Calculations and drawings. If the bounds is lower, use the second one..
ft=[a, b, c, d, e, f, g];
HighB1=1+max(abs(ft(2:end)))/abs(ft(1));
% Looks like the second one is the first inverted and with minus on the odd
% parts
ft2=[g, -f, e, -d, c, -b, a];
HighB2=1+max(abs(ft2(2:end)))/abs(ft2(1));
SndUsed = 0;
if HighB2 < HighB1
    ft=ft2;
    SndUsed = 1;
end
tic
t=SturmRoot(ft);
toc
if isempty(t)
    disp('No solutions');
    return;
end

x=linspace(-4,4);
y=polyval(ft,x);
plot(x,y);
y=polyval(ft2,x);
plot(x,y);
% Result will be f(theta). For each theta found, will find x,y
for i=1:numel(t)
    % For each then find B2 and B3
    theta=2*atan(t(i));
    if SndUsed==1
        theta = theta + pi;
    end
    thetaDeg = theta*180/pi
    a2=side*cos(theta)-1;
    b2=side*sin(theta);
    a3=side*cos(pi/3+theta);
    b3=side*sin(pi/3+theta)-1;
    den=2*(a2*b3-a3*b2); % assert != 0
    x=(b3*(L2^2-L1^2-a2^2-b2^2)-b2*(L3^2-L1^2-a3^2-b3^2))/den;
    y=(-a3*(L2^2-L1^2-a2^2-b2^2)+a2*(L3^2-L1^2-a3^2-b3^2))/den;
    
    B1=[x;y];
    B2=[x+side*cos(theta);y+side*sin(theta)];
    B3=[x+side*cos(pi/3+theta); y+side*sin(pi/3+theta)];
    figure
    fill( [B1(1) B2(1) B3(1)], [B1(2) B2(2) B3(2)], 'r' );
    hold on;
    plot([A1(1) B1(1)], [A1(2) B1(2)]);
    plot([A2(1) B2(1)], [A2(2) B2(2)]);
    plot([A3(1) B3(1)], [A3(2) B3(2)]);
    axis equal;
end

end