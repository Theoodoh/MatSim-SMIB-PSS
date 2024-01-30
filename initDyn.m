%% This m-file initiats MatSim1.0
clc
clear all
close all

%% Define Variables
global_var

%% User Must Set This path based on his/her PC
addpath('C:\Users\EBUKA\Desktop\matpower4.1')

%% Execute Power Flow Using MatPower4.1
pf = runpf(caseSMIB);

Vm=pf.bus(:,8);
Va=pf.bus(:,9);
Vb=Vm.*exp(1j*(Va*pi/180));

sys_data = caseSMIB;
baseMVA=sys_data.baseMVA;
bus=sys_data.bus;
branch1=sys_data.branch;
BRX=branch1(:,1:5);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch1);
YB=full(Ybus);
clear Ybus Yf Yt

v=Vb;

Pg=pf.gen(:,2)/baseMVA;
Qg=pf.gen(:,3)/baseMVA;
Pd=sys_data.bus(:,3)/baseMVA;
Qd=sys_data.bus(:,4)/baseMVA;
Ybus=YB;

%% Enter Dynamic parameters for the Generators
m=1;
n=2;

%         H      Xd     Xpd      Xq     Xpq    Tpd0    Tpq0
GenData=[ 6.4   0.8958 0.1198  0.8645  0.1969  6      0.535];


%% Enter Dynamic parameters for the Excitation Systems
%        KA  TA   KE   TE     KF       TF
ExData=[ 50  0.05  1    0.314  0.063    0.35];

%%
H=GenData(:,1);
Xd=GenData(:,2);
Xpd=GenData(:,3);
Xq=GenData(:,4);
Xpq=GenData(:,5);
Tpd0=GenData(:,6);
Tpq0=GenData(:,7);
ws=377;
M=2*H/ws;
D(1,1)=0.1*M(1);
D=0.2*M;
D=0.0;
KA=ExData(:,1);
TA=ExData(:,2);
KE=ExData(:,3);
TE=ExData(:,4);
KF=ExData(:,5);
TF=ExData(:,6);
Rs=[0];

YN=Ybus;

%% apply fault on the system
% fault=1 ; the system with fault
% fault=0 ; the system without fault

fault=1;

%% Select fault Type
% ftype=1; small change in references
% ftype=2; symetrical three phase faults
% ftype=3; load changes

ftype=2;


fPd=Pd;
fQd=Qd;

if ftype == 3
    fPd(5,1)=Pd(5,1)/2;
    fQd(5,1)=Qd(5,1)/2;
end
event=[0 1 1.1 10;
    0 0 0   0];
small_d=[
    0  1 1.1 10
    1  1.00 1 1;];
faulut_3ph=[
    0  1 1.1 10
    1  1.00 1 1;];
fYN = YN;
for i=1:n
    YL(i,i) = (Pd(i)-sqrt(-1)*Qd(i))/abs(v(i))^2;
    fYL(i,i) = (fPd(i)-sqrt(-1)*fQd(i))/abs(v(i))^2;
    
    YN(i,i) = YB(i,i)+YL(i,i);
    fYN(i,i) = YB(i,i)+fYL(i,i);
end
% YGG = YN(1:m,1:m);
% YGL = YN(1:m,m+1:n);
% YLG = YN(m+1:n,1:m);
% YLL = YN(m+1:n,m+1:n);
%
% Yred = YGG-YGL*inv(YLL)*YLG;
% Zred = inv(Yred);
% iYLL = inv(YLL);
Yred = YN;
if fault
    if ftype == 1   % small change in references
        fYred=Yred;
        small_d=[
            0  1 1.1 10
            1  1.05 1 1;];
    end
    if ftype == 2 % symetrical three phase faults
        fYred=Yred;
        faulut_3ph=[
            0  1 1.1 10
            1  0.00 1 1;];
    end
    if ftype == 3 % load changes
        
        fYred = Yred;
        
    end
    if ftype == 4 % line outage
        
    end
end


% YLV=-inv(YN(m+1:n,m+1:n))*YN(m+1:n,1:m);
%Initial Condition Calculation----------------------------------------
V0=abs(v);
Teta0=angle(v);

% for i=1:m
IG0=(Pg(2)-1j*Qg(2))/conj(v(2));
delta0=angle(v(2)+(Rs+1j*Xq)*IG0);

Id0=real(abs(IG0)*exp(1j*(pi/2+angle(IG0)-delta0)));
Iq0=imag(abs(IG0)*exp(1j*(pi/2+angle(IG0)-delta0)));
Vd0=real(V0(2)*exp(1j*(pi/2+Teta0(2)-delta0)));
Vq0=imag(V0(2)*exp(1j*(pi/2+Teta0(2)-delta0)));

Epq0=Vq0+Rs*Iq0+Xpd*Id0;
Epd0=(Xq-Xpq)*Iq0;
Efd0=Epq0+(Xd-Xpd)*Id0;
TM=Epq0*Iq0+(Xq-Xpd)*Id0*Iq0;
VR0=(KE+0.0039*exp(1.555*Efd0))*Efd0;
Rf0=KF*Efd0/TF;
Vref=V0(2)+VR0/KA;
Vref1=V0(2)+Efd0/KA;
% end
%%   ALgebraic equation

Rsm=diag(Rs);
Xpqm=diag(Xpq);
Xpdm=diag(Xpd);

XXpdq=[Xpq-Xpd];
XXdpd=[Xd-Xpd];
XXqpq=[Xq-Xpq];

Zdq=[Rsm -Xpqm;Xpdm Rsm];
iZdq=inv(Zdq);
Vinf =V0(1);
% Yred*v
% 
% Xs = 0.0625;
% alfa0 = asin(Pg(2)*Xs/V0(1)/V0(2));
% 
% IG = (V0(2)*exp(1j*alfa0)-V0(1))/(1j*Xs)