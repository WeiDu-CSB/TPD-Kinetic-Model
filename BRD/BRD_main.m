%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import BRD time course data reconstructed from Riching et al. fit 
BRD_tc=importdata('BRD_tc.csv');
BRD_data=BRD_tc.data(2:43,:);
BRD_data_time=BRD_tc.data(1,:);
BRD_data_anno=BRD_tc.textdata(2:43,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

% constants
V=4*10^(-12);   % cell volumn; L
NA_MOLE=6.02*10^23;  % Avocado's number

% target protein and E3 copy number; typical range: 1000 - 100000
N_CRBN=12265;   % CRBN protein copy number
N_VHL=15289;    % VHL protein copy number
N_BRD2=12162;   % BRD2 protein copy number
N_BRD3=37040;   % BRD3 protein copy number
N_BRD4=27613;   % BRD4 protein copy number

% target half-life; typical range: 0.5-200 hr
Thalf_BRD2=20.8;   % BRD2 protein half-life; hr
Thalf_BRD3=52.2;   % BRD3 protein half-life; hr
Thalf_BRD4=41;     % BRD4 protein half-life; hr

% degrader-target KD; typical range: 0.001-1 uM
KD_BRD2=0.228;  % dBET1/MZ1 BRD2 KD; uM
KD_BRD3=0.115;  % dBET1/MZ1 BRD3 KD; uM
KD_BRD4=0.120;  % dBET1/MZ1 BRD4 KD; uM

% degrader-E3 KD; typical range: 0.001-10 uM 
KD_CRBN=1.795;  % dBET1 CRBN KD; uM
KD_VHL=0.347;   % MZ1 VHL KD; uM

% cooperativity; typical range: 0.1-100
alpha_dBET1_BRD2=0.2;  % dBET1 BRD2 alpha; 
alpha_dBET1_BRD3=0.2;  % dBET1 BRD3 alpha; 
alpha_dBET1_BRD4=0.2;  % dBET1 BRD4 alpha; 
alpha_MZ1_BRD2=2.9;    % MZ1 BRD2 alpha;
alpha_MZ1_BRD3=11;     % MZ1 BRD3 alpha;
alpha_MZ1_BRD4=18;     % MZ1 BRD4 alpha;

% degrader competitive displacement IC50 for intracellular concentration
% calculation; typcial range: 0.01-1 uM
IC50_JQ1=0.054;     % JQ1 competitive displacement IC50; uM
IC50_dBET1=0.747;   % dBET1 competitive displacement IC50; uM
IC50_MZ1=0.432;     % MZ1 competitive displacement IC50; uM

% timespan of degradation; typical range: 30-2160 min
T=5*60;              % timespan of degradation; min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter conversion
T_CRBN=N_CRBN/V/NA_MOLE*10^6;   % CRBN protein concentration; uM
T_VHL=N_VHL/V/NA_MOLE*10^6;     % VHL protein concentration; uM
T_BRD2=N_BRD2/V/NA_MOLE*10^6;   % BRD2 protein concentration; uM
T_BRD3=N_BRD3/V/NA_MOLE*10^6;   % BRD3 protein concentration; uM
T_BRD4=N_BRD4/V/NA_MOLE*10^6;   % BRD4 protein concentration; uM
  
kdp_BRD2=log(2)/Thalf_BRD2/60;   % BRD2 endogeneous protein degradation rate; min^(-1)
kdp_BRD3=log(2)/Thalf_BRD3/60;   % BRD3 endogeneous protein degradation rate; min^(-1)
kdp_BRD4=log(2)/Thalf_BRD4/60;   % BRD4 endogeneous protein degradation rate; min^(-1)

ksp_BRD2=kdp_BRD2*T_BRD2;   % BRD2 endogeneous protein synthesis rate; uM/min;      
ksp_BRD3=kdp_BRD3*T_BRD3;   % BRD3 endogeneous protein synthesis rate; uM/min; 
ksp_BRD4=kdp_BRD4*T_BRD4;   % BRD4 endogeneous protein synthesis rate; uM/min; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intracellular dose calculation
dose_ori=[0.001 0.004 0.012 0.037 0.111 0.333 1];  % extracellular administrated concentration; uM
dose_dBET1=10.^(log10(IC50_JQ1*1000)/log10(IC50_dBET1*1000)*log10(dose_ori*1000))/1000;  % dBET1 intracellular concentration; uM
dose_dJQ1=10.^(log10(IC50_JQ1*1000)/log10(IC50_MZ1*1000)*log10(dose_ori*1000))/1000;  % JQ1 intracellular concentration; uM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degradation rate fitting
T_P=[T_BRD2 T_BRD3 T_BRD4];
T_E3=[T_CRBN T_VHL];
KD_P=[KD_BRD2 KD_BRD3 KD_BRD4];
KD_E3=[KD_CRBN KD_VHL];
alpha=[alpha_dBET1_BRD2 alpha_dBET1_BRD3 alpha_dBET1_BRD4 alpha_MZ1_BRD2 alpha_MZ1_BRD3 alpha_MZ1_BRD4];
k_dp=[kdp_BRD2 kdp_BRD3 kdp_BRD4];
k_sp=[ksp_BRD2 ksp_BRD3 ksp_BRD4];
dose=[dose_dBET1; dose_dJQ1];

DMAX_data=reshape(BRD_data(:,501),7,6)';   % initialize DMAX data matrix

kpr_series=10.^[-2:0.01:2];                % define kpr series 
residual=zeros(6,length(kpr_series));      % initialize residual

options=odeset('events',@events_time);    % set ODE timepoints
tic;

% simulate degradation time course at various kpr values
for k=1:length(kpr_series)

    disp(k)
    DMAX_fit=zeros(6,7);

    for i=1:2
  
       for j=1:7

          par=[KD_P T_E3(i) KD_E3(i) alpha(3*i-2) alpha(3*i-1) alpha(3*i) k_dp k_sp kpr_series(k) kpr_series(k) kpr_series(k) dose(i,j)];
          init=T_P;

          [t,y] = ode45(@BRD_model,[0 T],init,options,par);

          DMAX_fit(3*i-2,j)=y(length(t),1)/init(1);
          DMAX_fit(3*i-1,j)=y(length(t),2)/init(2);
          DMAX_fit(3*i,j)=y(length(t),3)/init(3);

       end

     end

     residual(:,k)= sum((DMAX_data-DMAX_fit).^2,2);

end

[min_value min_index]=min(residual');     % select best fit kpr value
k_pr=kpr_series(min_index);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time course (Figure 2A)

linecolor_1=bone(8);
linecolor_2=pink(8);
linecolor=[linecolor_1(flip(1:7),:);linecolor_2(flip(1:7),:)];
legend_cord=[[0.75 0.7 0.08 0.08];[0.75 0.25 0.08 0.08]];
annotation_cord=[[0.768 0.87 0 0];[0.775 0.42 0 0]];
protac=["dBET1" "MZ1"];
DMAX=zeros(6,7);
figure;
set(gcf,'position',[1200 400 1200 400]);

for i=1:2 
  
    for j=1:7
     
     par=[KD_P T_E3(i) KD_E3(i) alpha(3*i-2) alpha(3*i-1) alpha(3*i) k_dp k_sp k_pr(3*i-2) k_pr(3*i-1) k_pr(3*i) dose(i,j)];
     init=T_P;

     [t,y] = ode45(@BRD_model,[0 T],init,options,par);

     subplot(2,4,4*i-3);
     plot(t/60,y(:,1)/init(1),'--','LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
     xlim([0 T/60]);
     ylim([0 1]);
     title('BRD2',FontSize=15);
     xlabel('Time (hr)');
     ylabel('Target Level');
     hold on;
     subplot(2,4,4*i-2);
     plot(t/60,y(:,2)/init(2),'--','LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
     xlim([0 T/60]);
     ylim([0 1]);
     title('BRD3',FontSize=15);
     xlabel('Time (hr)');
     ylabel('Target Level');
     hold on;
     subplot(2,4,4*i-1);
     plot(t/60,y(:,3)/init(3),'--','LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
     xlim([0 T/60]);
     ylim([0 1]);
     title('BRD4',FontSize=15);
     xlabel('Time (hr)');
     ylabel('Target Level');
     hold on;

     DMAX(3*i-2,j)=y(length(t),1)/init(1);
     DMAX(3*i-1,j)=y(length(t),2)/init(2);
     DMAX(3*i,j)=y(length(t),3)/init(3);

    end

    subplot(2,4,4*i-3);
    for j=1:7
       p1=plot(BRD_data_time,BRD_data(i*21-20+j-1,:),'LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
       xlim([0 T/60]);
       ylim([0 1]);
    end
    hold off;
    subplot(2,4,4*i-2);
    for j=1:7
      p2=plot(BRD_data_time,BRD_data(i*21-13+j-1,:),'LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
      xlim([0 T/60]);
      ylim([0 1]);
    end
    hold off;
    subplot(2,4,4*i-1);
    for j=1:7
      p3=plot(BRD_data_time,BRD_data(i*21-6+j-1,:),'LineWidth',2.5,'Color',linecolor(7*i-7+j,:));
      xlim([0 T/60]);
      ylim([0 1]);
    end
    hold off;

  lgd=legend('0.001 uM','0.004 uM','0.012 uM','0.037 uM','0.111 uM','0.333 uM','1 uM','0.001 uM','0.004 uM','0.012 uM','0.037 uM','0.111 uM','0.333 uM','1 uM','Location','bestoutside','NumColumns',2);
  lgd.Position=legend_cord(i,:);
  title(lgd,[ "Model             Data"]);
  annotation('textbox',annotation_cord(i,:),'String',protac(i),'VerticalAlignment','Bottom','Edgecolor','none',FontSize=15);

end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
saveas(gcf,'Figure2A.pdf');


% plot DMAX (Figure 2B)

DMAX_color=lines(3);
legend_cord_2C=[[0.6 0.7 0.08 0.08];[0.6 0.25 0.08 0.08]];

figure;
set(gcf,'position',[600 500 600 500]);

for i=1:2

    subplot(2,2,i*2-1);

    for j=1:3
    
        semilogx(dose_ori,DMAX(3*i-2+j-1,:)*100,'+',MarkerSize=8,Color=DMAX_color(j,:),LineWidth=2);
        xlim([0.001 1]);
        ylim([0 100]);
        title(protac(i),FontSize=15);
        xlabel('Dose (uM)');
        ylabel('DMAX');
        ytickformat('percentage');
        hold on;
        semilogx(dose_ori,BRD_data((i*21-20+j*7-7):(i*21-20+j*7-1),501)*100,'o',MarkerSize=8,Color=DMAX_color(j,:),Linewidth=2);
        xlim([0.001 1]);
        ylim([0 100]);
        title(protac(i),FontSize=15);
        xlabel('Dose (uM)');
        ylabel('DMAX');
        ytickformat('percentage');
        hold on;

    end
    hold off;

    lgd=legend('BRD2','BRD2','BRD3','BRD3','BRD4','BRD4','Location','bestoutside','NumColumns',2,'Orientation','horizontal');
    lgd.Position=legend_cord_2C(i,:);
    title(lgd,[ "Model             Data"]);

end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
saveas(gcf,'Figure2B.pdf');



