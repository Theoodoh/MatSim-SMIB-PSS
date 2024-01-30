%% this m-file perfrmes small-signal analysis
initDyn

% AEO
% x = [39.6109 0.098433 0.0200001 0.0948553 0.0200001];
% x = [49.9999     0.813551    0.748044    0.197289    0.0200019];
% x = [5.4653 0.17051 0.14481 0.49772 0.02];
% x =[39.7285     0.118874         0.02    0.0792736    0.0200001];
% x = [30.6071    0.0776836     0.020563     0.135559    0.0224089];
x = [39.6109 0.098433 0.0200001 0.0948553 0.0200001];
KG = x(1);
Tw = 10;
T1 = x(2);
T2 = x(3);
T3 = x(4);
T4 = x(5);
Kpss = KG*T1*T3/(T2*T4);

%% Linearize Power System
% f11=linmod('SMIB');
f11=linmod('SMIB_pss');

% dx/dt = A.x + B.u
% y = C.x + D.u

Asys = f11.a ;
Bsys = f11.b ;
Csys = f11.c ;
Dsys = f11.d ;

%% Calculate Eigenvalues
egs = eig(Asys)
Ns = length(egs);

Damp = -real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq = abs(imag(egs))/(2*pi)

%% calculae Participation Factors
[Vs,D_eig] = eig(Asys);
Ws=inv(Vs);
for i=1:Ns
    for k=1:Ns
        Pfact1(k,i)=abs(Vs(k,i))*abs(Ws(i,k));
    end
end

for i=1:Ns
     Pfact(i,:)=Pfact1(i,:)/sum(Pfact1(i,:));
end

for i=1:Ns
    [s_val s_idx] = sort(Pfact(:,i),'descend');
    mod_idx(i,:) = s_idx(1:4)';
    pf_fact(i,:) = s_val(1:4)';
end
mod_idx;
pf_fact;
