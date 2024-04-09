function dy=BRD_model(t,y,par)

dy = zeros(3,1);

KD_BRD2=par(1);
KD_BRD3=par(2);
KD_BRD4=par(3);
T_E=par(4);
K_E=par(5);
alpha_BRD2=par(6);
alpha_BRD3=par(7);
alpha_BRD4=par(8);
kdp_BRD2=par(9);
kdp_BRD3=par(10);
kdp_BRD4=par(11);
ksp_BRD2=par(12);
ksp_BRD3=par(13);
ksp_BRD4=par(14);
kpr_BRD2=par(15);
kpr_BRD3=par(16);
kpr_BRD4=par(17);
L=par(18);

P_BRD2=y(1);
P_BRD3=y(2);
P_BRD4=y(3);

f=P_BRD2+T_E+1/(alpha_BRD2*L)*(L+KD_BRD2)*(L+K_E);
PLE_BRD2=(f-sqrt(f^2-4*P_BRD2*T_E))/2;
f=P_BRD3+T_E+1/(alpha_BRD3*L)*(L+KD_BRD3)*(L+K_E);
PLE_BRD3=(f-sqrt(f^2-4*P_BRD3*T_E))/2;
f=P_BRD4+T_E+1/(alpha_BRD4*L)*(L+KD_BRD4)*(L+K_E);
PLE_BRD4=(f-sqrt(f^2-4*P_BRD4*T_E))/2;

dy(1)=ksp_BRD2-kdp_BRD2*P_BRD2-kpr_BRD2*PLE_BRD2;
dy(2)=ksp_BRD3-kdp_BRD3*P_BRD3-kpr_BRD3*PLE_BRD3;
dy(3)=ksp_BRD4-kdp_BRD4*P_BRD4-kpr_BRD4*PLE_BRD4;
