function dy=WEE1_model(t,y,par)

dy = zeros(2,1);

T_P=par(1);
T_E=par(2);
K_E=par(3);
alpha=par(4);
k_dp=par(5);
k_pr=par(6);
K_P=par(7);
L=par(8);
T_S=par(9);
k_act=par(10);
k_inact=par(11);

k_sp=k_dp*T_P;

K_i=K_P;


P=y(1);
S=y(2);

f=P+T_E+1/(alpha*L)*(L+K_P)*(L+K_E);
PLE=(f-sqrt(f^2-4*P*T_E))/2;

dy(1)=k_sp-k_dp*P-k_pr*PLE;
dy(2)=1/(1+L/K_i)*k_act*P*(T_S-S)-k_inact*S; 
