%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

% constants
V=4*10^(-12);   % cell volumn; L
NA_MOLE=6.02*10^23;  % Avocado's number

% target protein and E3 copy number; typical range: 1000 - 100000
N_CRBN=12265;   % CRBN protein copy number
N_VHL=15289;    % VHL protein copy number

% degrader-E3 KD; typical range: 0.001-1 uM
KD_CRBN=1.795;  % PROTAC CRBN KD; uM
KD_VHL=0.347;   % PROTAC VHL KD; uM

% parameter conversion
T_CRBN=N_CRBN/V/NA_MOLE*10^6;   % CRBN protein concentration; uM
T_VHL=N_VHL/V/NA_MOLE*10^6;     % VHL protein concentration; uM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_E3=[T_VHL T_CRBN];
KD_E3=[KD_VHL KD_CRBN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model fitting (at 1uM) and prediction (at 100nM) of Foretinib CRBN and VHL PROTACs

% import Foretinib data including degradation fold changes, degrader-target KDs,target protein half-lives, target protein copy numbers 
degradation=importdata('Foretinib_data.csv');  
degradation_data=degradation.data(:,1:4);          % target fold change; 
affinity_data=degradation.data(:,5:6)/1000;        % target KD; uM
halflife=degradation.data(:,7);                    % target protein half-life; hr
copynumber=degradation.data(:,8);                  % target protein copy number; 
kdp=log(2)/halflife/60;                            % target endogeneous protein degradation rate; min^(-1)
concentration=copynumber/V/NA_MOLE*10^6;           % target protein concentration; uM

dose_1=1;     % PROTAC dose; uM
dose_2=0.1;   % PROTAC dose; uM

T=24*60;   % timespan; min

options=odeset('events',@events_time);       % set ODE timepoints
tic;

kpr_series=10.^[-2:0.01:2];       % define kpr series

min_kpr=zeros(length(halflife),2);      % initialize best fit kpr
min_fc=zeros(length(halflife),4);       % initialize best fit fold change

% simulate degradation fold change at various kpr values to find best fit
% kpr at 1uM then predict degradation fold change at 0.1uM
for j=1:2
    for i=1:length(halflife)

        disp([i,j])

        if(degradation_data(i,j*2)~=1)

           residual=zeros(1,length(kpr_series));
           fc=zeros(1,length(kpr_series));

           for k=1:length(kpr_series)

              par=[concentration(i) kdp(i) kpr_series(k) affinity_data(i,j) T_E3(j) KD_E3(j) dose_1]; 

              init=concentration(i);

              [t,y] = ode45(@kinome_model,[0 T],init,options,par);

              fc(k)=y(length(t))/init;
              residual(k)=(y(length(t))/init-degradation_data(i,j*2))^2;

           end

           [min_value min_index]=min(residual);
           min_kpr(i,j)=kpr_series(min_index);
           min_fc(i,2*j)=fc(min_index);

           par=[concentration(i) kdp(i) kpr_series(min_index) affinity_data(i,j)  T_E3(j) KD_E3(j) dose_2]; 

           init=concentration(i);

           [t,y] = ode45(@kinome_model,[0 T],init,options,par);

           min_fc(i,2*j-1)=y(length(t))/init;

        else
           min_kpr(i,j)=NaN;
           min_fc(i,2*j-1)=NaN;
           min_fc(i,2*j)=NaN;

        end

     end
end

% output predicted kpr and degradation fold change
csvwrite('Foretinib_min_kpr.csv',min_kpr);
csvwrite('Foretinib_min_fc.csv',min_fc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model fitting (at 1uM) and prediction (at 100nM) of TAE684 CRBN PROTACs

% import TAE684 data including degradation fold changes, degrader-target KDs,target protein half-lives, target protein copy numbers 
TAE684=importdata('TAE684_data.csv');
degradation_data=TAE684.data(:,2);         % degradation fold change
affinity_data=TAE684.data(:,3)/1000;       % target KD; uM
halflife=TAE684.data(:,4);                 % target protein half-life; hr
copynumber=TAE684.data(:,5);               % target protein copy number; 
kdp=log(2)/halflife/60;                    % target endogeneous protein degradation rate; min^(-1)
concentration=copynumber/V/NA_MOLE*10^6;   % target protein concentration; uM

dose_1=1;     % PROTAC dose; uM
dose_2=0.1;   % PROTAC dose; uM

T1=5*60;      % timespan; min
T2=4*60;      % timespan; min

options=odeset('events',@events_time);       % set ODE timepoints
tic;

kpr_series=10.^[-2:0.01:2];       % define kpr series

min_kpr=zeros(length(halflife),1);      % initialize best fit kpr
min_fc=zeros(length(halflife),2);       % initialize best fit fold change

% simulate degradation fold change at various kpr values to find best fit
% kpr at 1uM then predict degradation fold change at 0.1uM
for i=1:length(halflife)

    disp(i)

           residual=zeros(1,length(kpr_series));
           fc=zeros(1,length(kpr_series));

           for k=1:length(kpr_series)

              par=[concentration(i) kdp(i) kpr_series(k) affinity_data(i) T_E3(2) KD_E3(2) dose_1]; 

              init=concentration(i);

              [t,y] = ode45(@kinome_model,[0 T1],init,options,par);

              fc(k)=y(length(t))/init;
              residual(k)=(y(length(t))/init-degradation_data(i))^2;

           end

           [min_value min_index]=min(residual);
           min_kpr(i)=kpr_series(min_index);
           min_fc(i,1)=fc(min_index);

           par=[concentration(i) kdp(i) kpr_series(min_index) affinity_data(i) T_E3(2) KD_E3(2) dose_2 ]; 

           init=concentration(i);

           [t,y] = ode45(@kinome_model,[0 T2],init,options,par);

           min_fc(i,2)=y(length(t))/init;

end

% output predicted kpr and degradation fold change
csvwrite('TAE684_min_kpr.csv',min_kpr);
csvwrite('TAE684_min_fc.csv',min_fc);

