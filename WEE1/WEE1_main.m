%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import WEE1 data
WEE1=importdata('WEE1_data.csv');
WEE1_VHL_fc=WEE1.data(1:7,1)/100;
WEE1_CRBN_fc=WEE1.data(8:20,1)/100;
pCDK1_VHL_fc=WEE1.data(1:7,2)/100;
pCDK1_CRBN_fc=WEE1.data(8:20,2)/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

% constants
V=4*10^(-12);          % cell volumn; L
NA_MOLE=6.02*10^23;    % Avocado's number

% target protein, substrate, and E3 copy number; typical range: 1000 - 100000
N_CRBN=12265;          % E3(CRBN) protein copy number
N_VHL=15289;           % E3(VHL) protein copy number
N_WEE1=2503;           % target(WEE1) protein copy number
N_CDK1=1462484;        % substrate(CDK1) protein copy number

% target half-life; typical range: 0.5-200 hr
Thalf_WEE1=2.73;       % target(WEE1) protein half-life; hr

% degrader-E3 KD; typical range: 0.001-10 uM 
KD_CRBN=1.8;           % PROTAC E3(CRBN) KD; uM
KD_VHL=0.185;          % PROTAC E3(VHL) KD; uM

% degrader-target KD; typical range: 0.001-10 uM 
KD_WEE1=0.013;         % PROTAC target(WEE1) KD; uM

% cooperativity; typical range: 0.1-100
alpha=1;               % ternary complex cooperativity

% timespan of degradation; typical range: 30-2160 min; should match with experiment, or set relatively long to examine attainment of steady-state
T=24*60;               % timespan of degradation; min

% kinase catalytic parameters
k_act=1000;            % kinase target(WEE1) catalytic rate; uM^(-1)min^(-1)   
r=307;                 % kinase substrate(CDK1) inhibition constant; (see definition and calculation method in SI)

dose=1;                % degrader dosage; uM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter conversion
T_CRBN=N_CRBN/V/NA_MOLE*10^6;    % E3(CRBN) protein concentration; uM
T_VHL=N_VHL/V/NA_MOLE*10^6;      % E3(VHL) protein concentration; uM
T_P=N_WEE1/V/NA_MOLE*10^6;       % target(WEE1) protein concentration; uM
kdp=log(2)/Thalf_WEE1/60;        % target(WEE1) endogeneous protein degradation rate; min^(-1)
T_S=N_CDK1/V/NA_MOLE*10^6;       % substrate(CDK1) protein concentration; uM
k_inact=k_act*T_P/r;             % kinase substrate(CDK1) cellular inactivation rate; min^(-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degradation rate fitting
options=odeset('events',@events_time);     % set ODE timepoints
tic;

kpr_series=10.^[-2:0.01:2];       % define kpr series
dose_series=10.^[-4:0.2:1];       % define dose series

VHL_fc_fit=zeros(length(WEE1_VHL_fc),length(kpr_series));     % initialize VHL-based degrader fold change fit
for i=1:length(WEE1_VHL_fc)
    
    disp(i)

    for j=1:length(kpr_series)

           par=[T_P T_VHL KD_VHL alpha kdp kpr_series(j) KD_WEE1 dose 0 0 0]; 
           init=[T_P 0];
           [t,y] = ode45(@WEE1_model,[0 T],init,options,par);
           VHL_fc_fit(i,j)=y(length(t),1)/init(1);

    end
end
residual=(VHL_fc_fit-WEE1_VHL_fc).^2;          %find best fit kpr values
[min_value VHL_min_index]=min(residual');
VHL_kpr_fit=kpr_series(VHL_min_index);


CRBN_fc_fit=zeros(length(WEE1_CRBN_fc),length(kpr_series));     % initialize CRBN-based degrader fold change fit
for i=1:length(WEE1_CRBN_fc)
    
    disp(i)

    for j=1:length(kpr_series)
 
           par=[T_P T_CRBN KD_CRBN alpha kdp kpr_series(j) KD_WEE1 dose 0 0 0]; 
           init=[T_P 0];
           [t,y] = ode45(@WEE1_model,[0 T],init,options,par);
           CRBN_fc_fit(i,j)=y(length(t),1)/init(1);

    end
end
residual=(CRBN_fc_fit-WEE1_CRBN_fc).^2;         %find best fit kpr values
[min_value CRBN_min_index]=min(residual');
CRBN_kpr_fit=kpr_series(CRBN_min_index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degradation and functional inhibition prediction
VHL_fc=zeros(length(WEE1_VHL_fc),2);
for i=1:length(WEE1_VHL_fc)
    
    disp(i)   

    par=[T_P T_VHL KD_VHL alpha kdp VHL_kpr_fit(i) KD_WEE1 dose T_S k_act k_inact]; 
    init=[T_P k_act*T_P*T_S/(k_act*T_P+k_inact)];
    [t,y] = ode45(@WEE1_model,[0 T],init,options,par);
    VHL_fc(i,1)=y(length(t),1)/init(1);
    VHL_fc(i,2)=y(length(t),2)/init(2);

end

CRBN_fc=zeros(length(WEE1_CRBN_fc),2);
for i=1:length(WEE1_CRBN_fc)
    
    disp(i)
    par=[T_P T_CRBN KD_CRBN alpha kdp CRBN_kpr_fit(i) KD_WEE1 dose T_S k_act k_inact]; 
    init=[T_P k_act*T_P*T_S/(k_act*T_P+k_inact)];
    [t,y] = ode45(@WEE1_model,[0 T],init,options,par);
    CRBN_fc(i,1)=y(length(t),1)/init(1);
    CRBN_fc(i,2)=y(length(t),2)/init(2);
 
end

WEE1_kpr=[VHL_kpr_fit CRBN_kpr_fit]';
WEE1_out=[WEE1_kpr [VHL_fc;CRBN_fc]];
csvwrite('WEE1_out.csv',WEE1_out);         % save prediction output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate functional inhibition landscape

kpr_series=10.^[-2:0.01:2];       % define kpr series
KD_series=10.^[-4:0.01:1];        % define degrader-target KD series
dose_series=10.^[-4:0.2:1];       % define dose series

IC50=zeros(length(KD_series),length(kpr_series));      % initialize IC50 data matrix
IMAX=zeros(length(KD_series),length(kpr_series));      % initialize IMAX data matrix

for i=1:length(KD_series)
    for j=1:length(kpr_series)

        disp([i,j])

        fc=zeros(1,length(dose_series));

        for k=1:length(dose_series)

          par=[T_P T_CRBN KD_CRBN alpha kdp kpr_series(j) KD_series(i) dose_series(k) T_S k_act k_inact]; 
          init=[T_P k_act*T_P*T_S/(k_act*T_P+k_inact)];
          [t,y] = ode23(@WEE1_model,[0 T],init,options,par);
          fc(k)=y(length(t),2)/init(2);

        end

        [min_fc min_fc_index]=min(fc);       % find IMAX and associated dose
        IMAX(i,j)=1-min_fc;

        if(min(fc)>0.5)                      % If IMAX<50%, IC50=10uM
           IC50(i,j)=1;
        else
            conc=[-4:0.001:log10(dose_series(min_fc_index))];         % find DC50 and associated dose
            fit=spline(log10(dose_series(1:min_fc_index)),fc(1:min_fc_index),conc);
            [min_value min_conc_index]=min((1-fit-0.5).^2);   
            IC50(i,j)=conc(min_conc_index);
        end

    end
end

% save predicted IMAX and IC50 landscape
csvwrite('WEE1_IMAX.csv',IMAX);
csvwrite('WEE1_IC50.csv',IC50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot inhibition landscape with measured degrader data (Figure 5A)

figure;
set(gcf,'position',[720 230 720 230]);
marker_color=gray(8);

subplot(1,2,1);
imagesc(1-flip(IMAX));
xlabel('Degradation Rate (min^{-1})');
ylabel('Target KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('IMAX (%)',FontSize=15);
colorbar('Ticks',[0:0.1:1],'TickLabels',{'100%' '90%' '80%','70%','60%','50%','40%','30%','20%','10%','0%'});
hold on;
plot(CRBN_min_index,502-212,"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold off;

subplot(1,2,2);
imagesc(flip(IC50));
xlabel('Degradation Rate (min^{-1})');
ylabel('Targert KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('IC50 (nM)',FontSize=15);
colorbar('Ticks',[-4:1:1],'TickLabels',{'10^{-1}' '10^{0}' '10^{1}','10^{2}','10^{3}','10^4'});
hold on;
plot(CRBN_min_index,502-212,"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold off;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure5A.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot comparison between predicted and measured degradation and inhibition efficacy (Figure 5B)

figure;
set(gcf,'position',[520 230 520 230]);

subplot(1,2,1);
plot(WEE1_VHL_fc,VHL_fc(:,1),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0 1]);
ylim([0 1]);
xlabel('Observed WEE1 FC');
ylabel('Fitted WEE1 FC ');
hold on;
plot(WEE1_CRBN_fc,CRBN_fc(:,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0 1]);
ylim([0 1]);
xlabel('Observed WEE1 FC');
ylabel('Fitted WEE1 FC');
hold off;
legend('VHL','CRBN','Location','northoutside','NumColumns',2,'Orientation','horizontal');

subplot(1,2,2);
plot(pCDK1_VHL_fc,VHL_fc(:,2),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0 1]);
ylim([0 1]);
xlabel('Observed pCDK1 FC');
ylabel('Predicted pCDK1 FC ');
hold on;
plot(pCDK1_CRBN_fc,CRBN_fc(:,2),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0 1]);
ylim([0 1]);
xlabel('Observed pCDK1 FC');
ylabel('Predicted pCDK1 FC');
hold off;
legend('VHL','CRBN','Location','northoutside','NumColumns',2,'Orientation','horizontal');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure5B.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate degradation and functional inhibition hook effect landscape

kpr_series=10.^[-2:0.01:2];         % define kpr series
KD_series=10.^[-4:0.01:1];          % define degrader-target KD series
dose_series=10.^[-4:0.2:1];         % define dose series

DMAX=zeros(length(KD_series),length(kpr_series));        % initialize DMAX data matrix
DHOOK=zeros(length(KD_series),length(kpr_series));       % initialize DHOOK data matrix
IHOOK=zeros(length(KD_series),length(kpr_series));       % initialize IHOOK data matrix

for i=1:length(KD_series)

     for j=1:length(kpr_series)

         disp([i,j])

         D_fc=zeros(1,length(dose_series));
         I_fc=zeros(1,length(dose_series));
        
         for k=1:length(dose_series)

             par=[T_P T_CRBN KD_CRBN alpha kdp kpr_series(j) KD_series(i) dose_series(k)  T_S k_act k_inact]; 
             init=[T_P k_act*T_P*T_S/(k_act*T_P+k_inact)];
             [t,y] = ode23(@WEE1_model,[0 T],init,options,par);

             D_fc(k)=y(length(t),1)/init(1);
             I_fc(k)=y(length(t),2)/init(2);

         end

         DMAX(i,j)=1-min(D_fc);
         DHOOK(i,j)=(1-min(D_fc))-(1-D_fc(k));
         IHOOK(i,j)=(1-min(I_fc))-(1-I_fc(k));

    end
end

% save predicted DHOOK and IHOOK landscape
csvwrite('WEE1_DHOOK.csv',DHOOK);
csvwrite('WEE1_IHOOK.csv',IHOOK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot degradation hook effect (Figure 5E) and DMAX

figure;
set(gcf,'position',[720 230 720 230]);
marker_color=gray(8);

subplot(1,2,1);
imagesc(flip(DHOOK));
xlabel('Degradation Rate (min^{-1})');
ylabel('Target KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('DHOOK (%)',FontSize=15);
colorbar('Ticks',[0:0.1:1],'TickLabels',{'0%' '10%' '20%','30%','40%','50%','60%','70%','80%','90%','100%'});

subplot(1,2,2);
imagesc(1-flip(DMAX));
xlabel('Degradation Rate (min^{-1})');
ylabel('Target KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('DMAX (%)',FontSize=15);
colorbar('Ticks',[0:0.1:1],'TickLabels',{'100%' '90%' '80%','70%','60%','50%','40%','30%','20%','10%','0%'});

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure5E.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
