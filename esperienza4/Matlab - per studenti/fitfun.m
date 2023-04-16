function err = fitfun (BEST_PARS)

global ts1 ys1 gege  ts2 ys2 conc 

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

subplot (1,2,1); semilogx(ts1,yfit,'y',ts1,yMbCO,'r',ts1,yMb,'m',ts1,ytra1,'g',ts1,ys1,'*')
legend ('fit','Mb:CO','Mb','tra','dati')
subplot (1,2,2); semilogx(ts2,yfit2,'y',ts2,yMbCO2,'r',ts2,yMb2,'m',ts2,ytra2,'g',ts2,ys2,'*')
legend ('fit','Mb:CO','Mb','tra','dati')

err = norm(ys1-yfit)+norm(ys2-yfit2);
return