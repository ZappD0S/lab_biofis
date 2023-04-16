% Nonlinear curve fit with MINUIT algorithm.
                                                                                             
global ts1 ys1 ts2 ys2 conc 

conc=17e-6;
load T20.dat -ascii;
ts1=T20(1:length(T20),1);
ys1=T20(1:length(T20),2);
ys1=ys1*conc;
clear T20;

load T20atm01.dat -ascii;
ts2=T20atm01(1:length(T20atm01),1);
ys2=T20atm01(1:length(T20atm01),2);
ys2=ys2*conc;
clear T20atm01;
% grafico dei dati 
semilogx(ts1,ys1,'o', ts2,ys2,'o');
title('Input data')
    
%inizializzazione parametri
load par.txt -ascii;
START_GUESS=par;

[BEST_PARS,errs,chi2,errmatrix] = fminuit('fitfun',START_GUESS,'run_mode');

% visualizzazione finale fitting
clf;

[Tf1,Yf1] = ode15s(@funMb1,[ts1(1) ts1(length(ts1))],[0 0 0], 1e-15);
yMbCO =(BEST_PARS(6))+(interp1(Tf1(:,1),Yf1(:,1),ts1,'spline'));
yMb = interp1(Tf1(:,1),Yf1(:,2),ts1,'spline');  
ytra1= interp1(Tf1(:,1),Yf1(:,3),ts1,'spline'); 
yfit=yMbCO+yMb+ytra1;

[Tf2,Yf2] = ode15s(@funMb2,[ts2(1) ts2(length(ts2))],[0 0 0], 1e-15);
yMbCO2 =(BEST_PARS(10))+(interp1(Tf2(:,1),Yf2(:,1),ts2,'spline'));
yMb2 = interp1(Tf2(:,1),Yf2(:,2),ts1,'spline');  
ytra2= interp1(Tf2(:,1),Yf2(:,3),ts2,'spline'); 
yfit2=yMbCO2+yMb2+ytra2;
semilogx(ts1,yfit,'m',ts1,ys1,'m*', ts2,yfit2,'m',ts2,ys2,'m*')
legend ('fit','dati'), title ('finale')

k_1=BEST_PARS(1);
kout=BEST_PARS(2); 
kin=BEST_PARS(3); 
kc=BEST_PARS(4); 
k_c=BEST_PARS(5);
MbCO=BEST_PARS(6);
Mb=BEST_PARS(7);
tra1=BEST_PARS(8);
CO=BEST_PARS(9);
MbCO2=BEST_PARS(10);
M2b=BEST_PARS(11);
tra2=BEST_PARS(12);
CO2=BEST_PARS(13);

hold off