function dy = funMb(t,y)
global gege
BEST_PARS=gege;
dy = zeros(3,1);    % a column vector

%definisco le equazioni differenziali: 
%d[MbCO]/dt=dX/dt=y(1);
%d[Mb]/dt=dY/dt=y(2);
%d[tra1]/dt=dZ/dt=y(3);

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
Mb2=BEST_PARS(11);
tra2=BEST_PARS(12);
CO2=BEST_PARS(13);

%EQUAZIONI DIFFERENZIALI
dy(1)=[k_c*(tra1+y(3))+(-k_1-kc-kout)*(MbCO+y(1))+kin*(Mb+y(2))*(CO+y(2))];
dy(2)=[kout*(MbCO+y(1))-kin*(CO+y(2))*(Mb+y(2))];
dy(3)=[kc*(MbCO+y(1))-k_c*(tra1+y(3))];

return