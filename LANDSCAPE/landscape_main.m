%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

% constants
V=4*10^(-12);           % cell volumn; L
NA_MOLE=6.02*10^23;     % Avocado's number

% target protein and E3 copy number; typical range: 1000 - 100000
N_CRBN=12265;           % CRBN protein copy number
N_BTK=22692;            % BTK protein copy number

% target half-life; typical range: 0.5-200 hr
Thalf_BTK=71.9;         % BTK protein half-life; hr

% degrader-E3 KD; typical range: 0.001-10 uM 
KD_E3=3.84;             % CRBN KD; uM

% cooperativity; typical range: 0.1-100
alpha=1;                % ternary complex cooperativity

% timespan of degradation; typical range: 30-2160 min; should match with experiment, or set relatively long to examine attainment of steady-state
T=24*60;                % timespan of degradation; min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter conversion
T_E3=N_CRBN/V/NA_MOLE*10^6;     % E3(CRBN) protein concentration; uM
T_P=N_BTK/V/NA_MOLE*10^6;       % target(BTK) protein concentration; uM
kdp=log(2)/Thalf_BTK/60;        % target(BTK) endogeneous protein degradation rate; min^(-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate degradability landscape
options=odeset('events',@events_time);      % set ODE timepoints
tic;

kpr_series=10.^[-2:0.01:2];         % define kpr series
KD_series=10.^[-4:0.01:1];          % define degrader-target KD series
dose_series=10.^[-4:0.2:1];         % define dose series

DMAX=zeros(length(KD_series),length(kpr_series));     % initialize DMAX data matrix
DC50=zeros(length(KD_series),length(kpr_series));     % initialize DC50 data matrix

for i=1:length(KD_series)
    for j=1:length(kpr_series)

        disp([i,j])

        fc=zeros(1,length(dose_series));

        for k=1:length(dose_series)

           par=[T_P T_E3 KD_E3 alpha kdp kpr_series(j) KD_series(i) dose_series(k)]; 

           init=T_P;

           [t,y] = ode23(@landscape_model,[0 T],init,options,par);

           fc(k)=y(length(t))/init;

        end

        [min_fc min_fc_index]=min(fc);      % find DMAX and associated dose
        DMAX(i,j)=1-min_fc;

        if(min(fc)>0.5)                     % If DMAX<50%，DC50=10uM
            DC50(i,j)=1;
        else
            conc=[-4:0.001:log10(dose_series(min_fc_index))];                          % find DC50 and associated dose
            fit=spline(log10(dose_series(1:min_fc_index)),fc(1:min_fc_index),conc);     
            [min_value min_conc_index]=min((1-fit-0.5).^2);   
            DC50(i,j)=conc(min_conc_index);
  
        end

     end
end

% save predicted DMAX and DC50 landscape
csvwrite('BTK_DMAX.csv',DMAX);
csvwrite('BTK_DC50.csv',DC50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map BTK PROTACs onto degradability landscape

% import degrader data
BTK_data=importdata('BTK_data.csv');
BTK_KD=BTK_data.data(:,1)/1000;         % BTK PROTACs KD (uM)
Ramos_DC50=BTK_data.data(:,2)/1000;     % BTK PROTACs Ramos Cell DC50 (uM)
Ramos_DMAX=BTK_data.data(:,3)/100;      % BTK PROTACs Ramos Cell DMAX
THP1_DC50=BTK_data.data(:,4)/1000;      % BTK PROTACs THP-1 Cell DC50 (uM)
THP1_DMAX=BTK_data.data(:,5)/100;       % BTK PROTACs THP-1 Cell DMAX

BTK_map=zeros(7+6,4);
for i=1:7
    [~,KD_index]=min(abs(BTK_KD(i)-KD_series));
    BTK_map(i,1)=KD_index;     
    [~,kpr_index]=min(abs(Ramos_DC50(i)-10.^DC50(KD_index,:)));    % map kpr by measured DC50
    BTK_map(i,2)=kpr_index;
    BTK_map(i,3)=10^DC50(KD_index,kpr_index);                      % predict DC50 given mapped kpr
    BTK_map(i,4)=DMAX(KD_index,kpr_index);                         % predict DMAX given mapped kpr
end

for i=2:7
    [~,KD_index]=min(abs(BTK_KD(i)-KD_series));                
    BTK_map(7+i-1,1)=KD_index; 
    [~,kpr_index]=min(abs(THP1_DC50(i)-10.^DC50(KD_index,:)));     % map kpr by measured DC50
    BTK_map(7+i-1,2)=kpr_index;
    BTK_map(7+i-1,3)=10^DC50(KD_index,kpr_index);                  % predict DC50 given mapped kpr
    BTK_map(7+i-1,4)=DMAX(KD_index,kpr_index);                     % predict DMAX given mapped kpr
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot degradability landscape with measured degrader data (Figure 4A)
figure;
set(gcf,'position',[720 230 720 230]);
marker_color=gray(8);

subplot(1,2,1);
imagesc(1-flip(DMAX));
xlabel('Degradation Rate (min^{-1})');
ylabel('Target KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('DMAX (%)',FontSize=15);
colorbar('Ticks',[0:0.1:1],'TickLabels',{'100%' '90%' '80%','70%','60%','50%','40%','30%','20%','10%','0%'});
hold on;
%plot(BTK_out.data(1:7,2),502-BTK_out.data(1:7,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
plot(BTK_map(1:7,2),502-BTK_map(1:7,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold on;
plot(BTK_map(8:13,2),502-BTK_map(8:13,1),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
hold off;

subplot(1,2,2);
imagesc(flip(DC50));
xlabel('Degradation Rate (min^{-1})');
ylabel('Targert KD (nM)');
xticks([1 101 201 301 401]);
xticklabels({'10^{-2}' '10^{-1}' '10^{0}' '10^{1}' '10^{2}'});
yticks([1 101 201 301 401 501]);
yticklabels({'10^{4}' '10^{3}' '10^{2}' '10^{1}' '10^{0}' '10^{-1}'});
title('DC50 (nM)',FontSize=15);
colorbar('Ticks',[-4:1:1],'TickLabels',{'10^{-1}' '10^{0}' '10^{1}','10^{2}','10^{3}','10^4'});
hold on;
plot(BTK_map(1:7,2),502-BTK_map(1:7,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold on;
plot(BTK_map(8:13,2),502-BTK_map(8:13,1),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
hold off;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure4A.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot comparison between predicted and measured degrader efficacy (Figure 4B)

figure;
set(gcf,'position',[580 260 580 260]);

subplot(1,2,1);
loglog(Ramos_DC50*1000,BTK_map(1:7,3)*1000,"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.1 10000]);
ylim([0.1 10000]);
xticks([0.1 1 10 100 1000 10000]);
xlabel('Observed DC50 (nM)');
ylabel('Fitted DC50 (nM)');
hold on;
loglog(THP1_DC50(2:7)*1000,BTK_map(8:13,3)*1000,"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.1 10000]);
ylim([0.1 10000]);
xticks([0.1 1 10 100 1000 10000]);
xlabel('Observed DC50 (nM)');
ylabel('Fitted DC50 (nM)');
hold off;
legend('Ramos','THP-1','Location','northoutside','NumColumns',2,'Orientation','horizontal');

subplot(1,2,2);
plot(Ramos_DMAX,BTK_map(1:7,4),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.5 1]);
ylim([0.5 1]);
xticks([0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
yticks([0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
xlabel('Observed DMAX (%)');
ylabel('Predicted DMAX (%)');
hold on;
plot(THP1_DMAX(2:7),BTK_map(8:13,4),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.5 1]);
ylim([0.5 1]);
xticks([0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
yticks([0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
xlabel('Observed DMAX (%)');
ylabel('Predicted DMAX (%)');
hold off;
legend('Ramos','THP-1','Location','northoutside','NumColumns',2,'Orientation','horizontal');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure4B.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


