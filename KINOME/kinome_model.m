function dy=kinome_model(t,y,par)

dy = zeros(1,1);

T_P=par(1);
k_dp=par(2);
k_pr=par(3);
K_P=par(4);
T_E=par(5);
K_E=par(6);
L=par(7);

k_sp=k_dp*T_P;
alpha=1;

P=y(1);

f=P+T_E+1/(alpha*L)*(L+K_P)*(L+K_E);
PLE=(f-sqrt(f^2-4*P*T_E))/2;

dy(1)=k_sp-k_dp*P-k_pr*PLE;
