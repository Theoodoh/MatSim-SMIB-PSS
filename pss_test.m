% function J = pss_test
% clc
% clear all
% close all

% AEO
%x = [39.6109     0.098433    0.0200001    0.0948553    0.0200001];
% x = [49.9999     0.813551    0.748044    0.197289    0.0200019];
% x = [8.0927    0.10375    0.082337    0.36438    0.020015];
% x =[5.3476    0.15755    0.13608    0.51788    0.020044];

% x = [19.4017 0.182881 0.0351571 0.0767505 0.0205659];
load('sys_IO','f11')
As = f11.a;
Bs = f11.b;
Cs = f11.c;
Ds = f11.d;

Tw = 10;
KG = x(1);
T1 = x(2);
T2 = x(3);
T3 = x(4);
T4 = x(5);
Kpss = KG*T1*T3/(T2*T4);

b = [KG*T1*T3*Tw (KG*T1*Tw + KG*T3*Tw) KG*Tw 0];
a = [T2*T4*Tw  (T2*T4 + T2*Tw + T4*Tw) (T2 + T4 + Tw) 1];

[A,B,C,D]= tf2ss(b,a);

Af = (A);
Bf = (B);
Cf = (C);
Df = (D);

Asys_1 = As + Bs*Df*Cs;
Asys_2 = Bs*Cf;
Asys_3 = Bf*Cs;
Asys_4 = Af + Bf*Ds*Cf;
Asys = [Asys_1 Asys_2;
    Asys_3 Asys_4];

egs = eig(Asys);

[z_val z_idx]=sort(abs(egs),'descend');
egs_new=egs;
egs_new(z_idx(end-1:end))=[];

%% unstable modes
ss_idx = find(real(egs_new)>0);
uss = egs_new(ss_idx);

%% EM modes
% Damp=-real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq = abs(imag(egs_new))/(2*pi);
em_idx = find(freq>0 & freq<3);

objf = max(real(egs_new(em_idx)))+sum(egs_new(ss_idx));

if isempty(objf)
    objf = max(real(egs_new));
end
J = objf;