clear all
close all
clc
%% DANE
n=300;                                          %[rpm]

a=0.250;                                        %[m]
b=0.150;                                        %[m]
c=0.050;                                        %[m]

dlugosc=0:0.001:(a+b+c);                        %[m]        %beam leanght

rB=0.120;                                       %[m]
rC=0.120;                                       %[m]

F1x=5000;                                       %[N] is F1
F1y=3000;                                       %[N]

F2x=5000;                                       %[N]

Msb=F1x*rB;                                     %[kNm]      
Msc=Msb;                                        %[kNm]      
F2y=Msc/rC;                                     %[N]        

%% force of reaction

%calculate a sysytem of eq
syms RAx RAy RDx RDy
row1=F2x-F1x+RAx+RDx==0;                        %X axs
row2=F2y-F1y+RAy+RDy==0;                        %Yb axs
% 0=RDz;                                        %Z axs
row3=(a+b+c)*RAy-(b+c)*F1y+c*F2y==0;            
row4=-(a+b+c)*RAx+(b+c)*F1x-c*F2x==0;          

sol = solve([row1, row2, row3, row4], [RAx, RAy, RDx, RDy]);
RAx = round(double(sol.RAx));
RAy = round(double(sol.RAy));
RDx = round(double(sol.RDx));
RDy = round(double(sol.RAy));
disp("Solution: ")
disp(sol)
%% wykresy sił tnących i momentów
% z-x
i = 1;
Mx(i) = 0;
My(i) = 0;
for x = 0.001:0.001:a+b+c
    i = i + 1;
    if(x <= a )
        Mx(i) = -RAx * x;
        My(i) = -RAy * x;
    elseif(x > a && x <= a+b)
        Mx(i) = -RAx * x + F1x * (x-a);
        My(i) = -RAy * x + F1y * (x-a);
    elseif(x > a+b && x <= a+b+c)
        Mx(i) = -RAx * x + F1x * (x-a) - F2x * (x-a-b);
        My(i) = -RAy * x + F1y * (x-a) - F2y * (x-a-b);
    else
    end
end

i = 1;
Tx(i) = 0;
Ty(i) = 0;
for x = 0.001:0.001:a+b+c
    i = i + 1;
    if(x <= a )
        Tx(i) = -RAx;
        Ty(i) = -RAy;
    elseif(x > a && x <= a+b)
        Tx(i) = -RAx + F1x;
        Ty(i) = -RAy + F1y;
    elseif(x > a+b && x <= a+b+c)
        Tx(i) = -RAx + F1x - F2x;
        Ty(i) = -RAy + F1y - F2y;
    else
        Tx(i) = -RAx + F1x - F2x - RDx;
        Ty(i) = -RAy + F1y - F2y - RDy;
    end
end

figure('Name','Shear Forces','NumberTitle','off');

subplot(2,1,1)
plot(dlugosc,Ty)
title("Z-Y Plane")
ylabel("Shear force [N]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(Ty)-200 max(Ty)+200])
grid on

subplot(2,1,2)
plot(dlugosc,Tx)
title("Z-X Plane")
ylabel("Shear force [N]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(Tx)-200 max(Tx)+200])
grid on



figure('Name','Bending Moments','NumberTitle','off');

subplot(2,1,1)
plot(dlugosc,My)
title("Z-Y Plane")
ylabel("Bending moment [Nm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(My)-200 max(My)+200])
grid on

subplot(2,1,2)
plot(dlugosc,Mx)
title("Z-X Plane")
ylabel("Bending moment [Nm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(Mx)-200 max(Mx)+200])
grid on


%% Sum of bending moments:

for x=1:length(dlugosc)
    Mg(x)= sqrt(Mx(x)^2 + My(x)^2);
end
figure('Name','Net Moments','NumberTitle','off');
plot(dlugosc,Mg)
title('Net Moments')
ylabel("Net moment [Nm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(Mg)-200 max(Mg)+200])
grid on

%% graph - torque
MSkrec = zeros(1,451);
MSkrec(251:401) = Msc;

figure('Name','Torsional moments','NumberTitle','off');
plot(dlugosc, MSkrec)
title('Torsional moments')
ylabel("Torsional moment [Nm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(MSkrec)-200 max(MSkrec)+200])
grid on


%% Stress hypotheses - reduced moment H-M-H

for x=1:length(shaftlength)
    Mred(x) = sqrt(Mbnet(x)^2 + 3/4 * Mt(x)^2);
end

figure('Name','Reduced moments','NumberTitle','off');
plot(dlugosc,Mred)
title('Reduced moment')
ylabel("Reduced moment [Nm]")
xlabel("x [m]")
xlim([min(shaftlength) max(dlugosc)])
ylim([min(Mred)-200 max(Mred)+200])
grid on

%% minimum cross section
kg = 75*10^6;
for x=1:length(dlugosc)
    dmin(x) =  nthroot((32*Mred(x))/(pi*kg),3);
end

figure('Name','Minimum Diameter','NumberTitle','off');
plot(dlugosc, dmin)
title('Minimum diameter')
ylabel("Minimum diameter [m]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
ylim([min(dmin)-0.001 max(dmin)+0.001])
grid on
%% key length
lkb = 2*Msb/(50*1*5*150);
lkc = 2*Msc/(50*1*5*150);
%% bearing
RA=sqrt(RAx^2 + RAy^2);     
RD=sqrt(RDx^2 + RDy^2);    
ft=1;
fd=1.4;
PA=RA*ft*fd;
PD=RD*ft*fd; 
L=10000;
epsilon=3;
bearings.CA=(nthroot(((L*60*n*(PA^3))/10^6),3))/1000;
bearings.CD=(nthroot(((L*60*n*(PD^3))/10^6),3))/1000;
disp(bearings)


d1 = 0.0396;
d2 = 0.0454;
d3 = 0.0419;
d4 = 0.0278;
sf = 1.2; % safety factor
g = 0.0120;
Mt = 600;
% gear diameters

Db = 0.05;
Lp = 1.4*Db;
g = 0.275*Db;
Dw = Db + g;

% key calculations
% ISO 2491A
Bb = 0.014;
Hb = 0.006;
% ISO 2491A
Bc = 0.014;
Hc = 0.006;

% steel C45
Msb = 600;
Msc = 600;
Re = 335000000;
Xe = 2;
Kr = Re/Xe;
Zrc = 0.7*Kr;
Fb = 2*Msb/Db;
Fc = 2*Msc/Db;
lb = 2*Fb/(Hb*Zrc);
lc = 2*Fc/(Hc*Zrc);

% FOR POINT A
% from SKF catalogue
% 7209 BE
% C = 35.8kN
% d = 45mm
% D = 85mm
% B = 19mm
% limiting speed 7500rpm

% FOR POINT D
% from SKF catalogue
% 7207 BE
% C = 29.1kN
% d = 35mm
% D = 72mm
% B = 17mm
% limiting speed 9000rpm

Da = 46;
Db = 46;
Dc = 50;
Dd = 36;
%% Deflection
% points of interest
x0 = 0;
x1 = 0.010;
x2 = 0.240; % diameter change
x3 = 0.250; % b
x4 = 0.400;
x5 = 0.410; % diameter change
x6 = 0.450; % c
% integration constants
I1 = pi*0.045^4/64; % x0 - x2
I2 = pi*0.050^4/64; % x2 - x5
I3 = pi*0.036^4/64; % x5 - x6
I = I1+I2+I3;
E = 210000000000;
C = 17.5926;
i = 1;
wy(i) = 0;
for x = 0.001:0.001:a+b+c
    i = i + 1;
    if(x <= x1 )
        
    elseif(x > x1 && x <= x2)
        
    elseif(x > x2 && x <= x3)
        
    elseif(x > x3 && x <= x4)
        
    elseif(x > x4 && x <= x5)
        
    elseif(x > x5 && x <= x6)
        
    else
        wy(i) = 0;
    end
end
wy = wy*1000;
figure('Name','Deflections','NumberTitle','off');
subplot(2,1,1)
plot(dlugosc,wy)
title('Deflections Y-Z')
ylabel("Deflections Y-Z [mm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
grid on

C = 41.6667;
i = 1;
wx(i) = 0;
for x = 0.001:0.001:a+b+c
    i = i + 1;
    if(x <= x1 )
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 + C*x);
    elseif(x > x1 && x <= x2)
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 + C*x);
    elseif(x > x2 && x <= x3)
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 + C*x);
    elseif(x > x3 && x <= x4)
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 +1/6* F1x * (x-x3)^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 +1/6* F1y * (x-x3)^3 + C*x);
    elseif(x > x4 && x <= x5)
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 +1/6* F1x * (x-x3)^3 - 1/6*F2x *(x-x4)^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 +1/6* F1y * (x-x3)^3 - 1/6*F2y *(x-x4)^3 + C*x);
    elseif(x > x5 && x <= x6)
        wx(i) = 1/(E*I)*(-1/6*RAx * x^3 +1/6* F1x * (x-x3)^3 - 1/6*F2x *(x-x4)^3 + C*x);
        wy(i) = 1/(E*I)*(-1/6*RAy * x^3 +1/6* F1y * (x-x3)^3 - 1/6*F2y *(x-x4)^3 + C*x);
    else
        wx(i) = 0;
        wy(i) = 0;
    end
end
wx = wx*1000;
subplot(2,1,2)
plot(dlugosc,wx)
title('Deflections X-Z')
ylabel("Deflections X-Z [mm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
grid on

wnet = sqrt(wx.^2+wy.^2);
figure('Name','Net deflections','NumberTitle','off');
plot(dlugosc,wnet)
title('Net Deflections')
ylabel("Net Deflections [mm]")
xlabel("x [m]")
xlim([min(dlugosc) max(dlugosc)])
grid on

%% critical speed
E=210*10^9;
C=3.5;
d1=0.045;
d2=0.050;
d3=0.036;
%moments of inertia
I1 = pi*d1^4/64; % x0 - x2
I2 = pi*d2^4/64; % x2 - x5
I3 = pi*d3^4/64; % x5 - x6
%density
p=7800;
%crossectional areas
A1=pi*(d1/2)^2;
A2=pi*(d2/2)^2;
A3=pi*(d3/2)^2;
%length of section
l1=0.240;
l2=0.170;
l3=0.040;

Ncrit=C/(2*pi)*sqrt(E*I1/(p*A1*l1^4)+E*I2/(p*A2*l2^4)+E*I3/(p*A3*l3^4));