global ts1 ys1 ts2 ys2 conc gege 

gege=BEST_PARS;

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

%Creo la matrice con le curve sperimentale, soluz numerica e fitting
z1=[ts1 ys1 yMbCO yMb ytra1 yfit];

%Creo la matrice con le curve sperimentale, soluz numerica e fitting
z1=[ts1 ys1 yMbCO yMb ytra1 yfit];
z2=[ts2 ys2 yMbCO2 yMb2 ytra2 yfit2];

%Creo la matrice con i residui
residui = norm(ys1-yfit)+norm(ys2-yfit2);

%Salvo le curve
save 'fit_T20.dat' 'z1' -ascii;
save 'fit_T20_01.dat' 'z2' -ascii;
save 'par_1e01.dat' 'BEST_PARS' -ascii;
save 'residui.dat' 'residui' -ascii;

