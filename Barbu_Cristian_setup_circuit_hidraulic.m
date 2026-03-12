%%
% Nume si prenume: Barbu Cristian-Emilian
%

clearvars
clc
close all;

%% Magic numbers (replace with received numbers)
m = 4; 
n = 8; 

%% Process data (fixed, do not modify)
c1 = (1000+n*300)/10000;
c2 = (1.15+2*(m+n/10)/20);
a1 = 2*c2*c1;
a2 = c1;
b0 = (1.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [2*(m/2+rand(1)*m/5); m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 10*c2/c1/1e4*1.5; % fundamental step size
Tfin = 30*c2/c1;

gain = 10;
umin = 0; umax = gain; % input saturation
ymin = 0; ymax = b0*gain/1.5; % output saturation

whtn_pow_in = 1e-6*5*(((m-1)*8+n/2)/5)/2*6/8; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(10); % input quantizer (DAC)

whtn_pow_out = 1e-5*5*(((m-1)*25+n/2)/5)*6/80*(0.5+0.3*(m-2)); % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = (m/2+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
u0 = 0;        % fixed
ust = 3;     % must be modified (saturation)
t1 = 10*c2/c1; % recommended 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("Barbu_Cristian2_circuit_hidraulic_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
title("Sistemul de baza");
shg

%% System identification
i4 = 13287;
i3 = 12470;
i2 = 6506;
i1 = 6460;

u0 = mean(u(i1:i2));
ust = mean(u(i3:i4));

y0 = mean(y(i1:i2));
yst = mean(y(i3:i4));

K = (yst-y0)/(ust-u0);

%%T1 (constanta dominanta)

i5 = 7121;
i6 = 9182;
t_regresie = t(i5:i6);
y_regresie = log(yst-y(i5:i6));
figure;
plot(t_regresie,y_regresie); %cu cat ma apropii de stationar -> zgomot mult
title("Grafic regresie liniara");

A_regresie = [sum(t_regresie.^2),sum(t_regresie);sum(t_regresie),length(t_regresie)];
B_regresie = [sum(y_regresie.*t_regresie); sum(y_regresie)];

theta = inv(A_regresie)*B_regresie

T1 = -1/theta(1)

i7 = 6992;
i8 = 7436;

Ti = t(i8)-t(i7)

T2vec = 0.1:0.1:3.19;

Fun = T1*T2vec.*log(T2vec)-T2vec*(Ti+T1*log(T1))+T1*Ti;
figure; 
plot(T2vec,Fun);
title("Grafic pentru determinarea lui T2");

T2 =1.5; %se realizeaza intersectia axei oY cu oX si cand y = 0, astfel se ia valoarea lui x pentru constanta T2
%% Validare
H = tf(K,[T1*T2,T1+T2,1])
ysim = lsim(H,u,t);%conditii initiale nule

figure;
plot(t,u,t,y,t,ysim);
title("Validare sistem->conditii initiale NULE");

%Spatiul starilor
A = [0 1 ; -1/T1/T2 , -(1/T1+1/T2)];
B = [0 ; K/T1/T2];
C = [1 0];
D = 0;

sys = ss(A,B,C,D);

ysim2 = lsim(sys,u,t,[y(1),2.41]); %cu conditii initiale nenule;

figure;
plot(t,u,t,y,t,ysim2);
title("Validare sistem->conditii initiale nenule");

J = 1/sqrt(length(t)*norm(y-ysim2));
eMPN = norm(y-ysim2)/norm(y-mean(y))*100; %calcul eroare medie patratica normata
