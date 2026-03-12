%%
% Nume si prenume: Barbu Cristia-Emilian
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 4;
n = 8;

%% Process data and experiment setup (fixed, do not modify)
Ts = 500e-6; % fundamental step size

u_star = 1.2+n*0.075;
delta = 0.125;
delta_spab = 0.075;

umin = -5; umax = 5; % input saturation
ymin = -100; ymax = 100; % output saturation

g = 9.81;
% pendulum parameters
M = 0.8-n/48;
l = 1.2-m/24;
b = 0.3+m/24;
% measurement
c1 = 180/pi;
c2 = 4+n/2;

% (theta0,omega0)
rng(m+10*n)
x0_slx = [(n+3)/50,(-1)^(n+1)*m/20];

% input white noise power and sampling time
whtn_pow_in = 1e-10*(Ts*1e4)/2; 
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(13); % input quantizer (DAC)

% output white noise power and sampling time
whtn_pow_out = 1e-3*Ts; 
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(13); % output quantizer (ADC)

meas_rep = round(7+n/2); % data acquisition hardware sampling limitation

%% Input setup (can be changed/replaced/deleted)
%Tfin = 50; % simulation duration
%t0 = Tfin/4;
%t1 = Tfin/2;

%timpul necesar aclimatizarii sistemului
t1 = 18;
tu = 1;
tpo = 1.7;

%calibrare sistem
N = 4;
p = round(tpo/N/Ts);

DeltaT = (2^N-1)*p*Ts*3; %perioada SPAB
[input_LUT_dSPACE,Tfin] = generate_input_signal(Ts,t1,DeltaT,N,p,u_star,delta,delta_spab);


%% Data acquisition (use t, u, y to perform system identification)
out = sim("pendul_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

subplot(211)
plot(t,u)
subplot(212)
plot(t,y)
shg

%% System identification

i1 = 53473;
i2 = 108872;
i3 = 118886;
i4 = 174900;

N = 11;

t_id = t(i1:N:i2);
u_id = u(i1:N:i2);
u_id = u_id - mean(u_id);%eliminare componente stationare DC
y_id = y(i1:N:i2);
y_id = y_id - mean(y_id);

t_vd = t(i3:N:i4);
u_vd = u(i3:N:i4);
u_vd = u_vd - mean(u_vd);
y_vd = y(i3:N:i4);
y_vd = y_vd - mean(y_vd);

figure;
%semnale identificare
subplot(2,2,1);
plot(t_id,u_id);
subplot(2,2,3);
plot(t_id,y_id);
%semnale validare
subplot(2,2,2);
plot(t_vd,u_vd);
subplot(2,2,4);
plot(t_vd,y_vd);


Ts = t_id(2)-t_id(1);
dat_id = iddata(y_id,u_id,Ts);
dat_vd = iddata(y_vd,u_vd,Ts);

%% IV4
%nA = 2;
%nB = 2;
%nd = 1;
model_iv4 = iv4(dat_id,[2,2,1]);
figure; resid(model_iv4,dat_vd);
figure; compare(model_iv4,dat_vd);
%dupa validare => iv4 nu e ok

%% ARX

model_arx = arx(dat_id,[2,2,1]);
figure; resid(model_arx,dat_vd);
figure; compare(model_arx,dat_vd);

%% ARMAX  - VALID 
%nC = 1,2,3,... luam valori pana iasa ok
model_armax = armax(dat_id,[2,2,2,1]);
figure; resid(model_armax,dat_vd);
figure; compare(model_armax,dat_vd);
% pentru nC = 2  se imbunatateste metoda
%% OE 

model_oe = oe(dat_id,[2,2,1]);
figure; resid(model_oe,dat_vd);
figure; compare(model_oe,dat_vd);

%% BJ - VALID 

nD_bj = 1;nB_bj = 2;nC_bj = 1;nF_bj = 2;nd_bj = 1;
model_bj = bj(dat_id,[nB_bj, nC_bj, nD_bj, nF_bj, nd_bj]);
figure; resid(model_bj,dat_vd);
figure; compare(model_bj,dat_vd);

%% SPATIUL STARILOR

model_n4sid = n4sid(dat_id,2)
figure; resid(model_n4sid,dat_vd);
figure; compare(model_n4sid,dat_vd);
%graficul la compare nu e bun =>  iau alt set de date etc;
% Nu mai iau, s-a validat cu SSEST!!!

%% SSEST - VALID 

model_ssest = ssest(dat_id,2)
figure; resid(model_ssest,dat_vd);
figure; compare(model_ssest,dat_vd);

%zpk(model_ssest)  => functia de transfer ne arata ca exista un zerou de
%faza minima -> zerou in semiplanul stang

 