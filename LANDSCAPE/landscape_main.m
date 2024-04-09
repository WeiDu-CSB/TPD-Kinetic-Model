% input parameters
N_CRBN=12265;         % E3(CRBN) protein copy number
N_BTK=22692;          % target(BTK) protein copy number
Thalf_BTK=71.9;       % target(BTK) protein half-life; hr
KD_E3=3.84;           % PROTAC E3(CRBN) KD; uM
alpha=1;              % ternary complex cooperativity
V=4*10^(-12);         % cell volumn; L
NA_MOLE=6.02*10^23;   % Avocado's number
T=24*60;              % timespan of degradation; min

% parameter conversion
T_E3=N_CRBN/V/NA_MOLE*10^6;     % E3(CRBN) protein concentration; uM
T_P=N_BTK/V/NA_MOLE*10^6;       % target(BTK) protein concentration; uM
kdp=log(2)/Thalf_BTK/60;        % target(BTK) endogeneous protein degradation rate; min^(-1)


% simulate degradability landscape

options=odeset('events',@events_time);
tic;

kpr_series=10.^[-2:0.01:2];
KD_series=10.^[-4:0.01:1];
dose_series=10.^[-4:0.2:1];

DMAX=zeros(length(KD_series),length(kpr_series));
DC50=zeros(length(KD_series),length(kpr_series));

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

        [min_fc min_fc_index]=min(fc);
        
        DMAX(i,j)=1-min_fc;

        conc=[-4:0.001:log10(dose_series(min_fc_index))];
        fit=spline(log10(dose_series(1:min_fc_index)),fc(1:min_fc_index),conc);
        [min_value min_conc_index]=min((1-fit-(1-min_fc)/2).^2);

        if(min(fc)>0.5)
            DC50(i,j)=1;
        else
            conc=[-4:0.001:log10(dose_series(min_fc_index))];
            fit=spline(log10(dose_series(1:min_fc_index)),fc(1:min_fc_index),conc);
            [min_value min_conc_index]=min((1-fit-0.5).^2);   

            DC50(i,j)=conc(min_conc_index);
  
        end

        end
end

csvwrite('BTK_DMAX.csv',DMAX);
csvwrite('BTK_DC50.csv',DC50);


% plot degradability landscape with measured degrader data (Figure 4A)

BTK_out=importdata('BTK_out.csv');
BTK_index=importdata('BTK_index.csv');

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
plot(BTK_index.data(1:7,2),502-BTK_index.data(1:7,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold on;
plot(BTK_index.data(8:13,2),502-BTK_index.data(8:13,1),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
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
plot(BTK_index.data(1:7,2),502-BTK_index.data(1:7,1),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
hold on;
plot(BTK_index.data(8:13,2),502-BTK_index.data(8:13,1),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
hold off;


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0.1 0.1 screenposition(3:4)+0.1],'PaperSize',[screenposition(3:4)+0.3]);
saveas(gcf,'Figure4A.pdf');



% plot comparison between predicted and measured degrader efficacy (Figure 4B)

figure;
set(gcf,'position',[580 260 580 260]);

subplot(1,2,1);
loglog(BTK_out.data(1:7,3),BTK_out.data(1:7,5),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.1 10000]);
ylim([0.1 10000]);
xticks([0.1 1 10 100 1000 10000]);
xlabel('Observed DC50 (nM)');
ylabel('Fitted DC50 (nM)');
hold on;
loglog(BTK_out.data(8:13,3),BTK_out.data(8:13,5),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.1 10000]);
ylim([0.1 10000]);
xticks([0.1 1 10 100 1000 10000]);
xlabel('Observed DC50 (nM)');
ylabel('Fitted DC50 (nM)');
hold off;
legend('Ramos','THP-1','Location','northoutside','NumColumns',2,'Orientation','horizontal');

subplot(1,2,2);
plot(BTK_out.data(1:7,4),BTK_out.data(1:7,6),"pentagram",MarkerFaceColor=marker_color(7,:),MarkerEdgeColor="black",MarkerSize=8);
xlim([0.5 1]);
ylim([0.5 1]);
xticks([0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
yticks([0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({'50%' '60%' '70%' '80%' '90%' '100%'});
xlabel('Observed DMAX (%)');
ylabel('Predicted DMAX (%)');
hold on;
plot(BTK_out.data(8:13,4),BTK_out.data(8:13,6),"pentagram",MarkerFaceColor=marker_color(3,:),MarkerEdgeColor="black",MarkerSize=8);
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


