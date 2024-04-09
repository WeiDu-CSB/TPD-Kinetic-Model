% input parameters

N_CRBN=12265;   % CRBN protein copy number
N_VHL=15289;    % VHL protein copy number

KD_CRBN=1.795;  % PROTAC CRBN KD; uM
KD_VHL=0.347;   % PROTAC VHL KD; uM

V=4*10^(-12);   % cell volumn; L
NA_MOLE=6.02*10^23;  % Avocado's number


% parameter conversion

T_CRBN=N_CRBN/V/NA_MOLE*10^6;   % CRBN protein concentration; uM
T_VHL=N_VHL/V/NA_MOLE*10^6;     % VHL protein concentration; uM

T_E3=[T_VHL T_CRBN];
KD_E3=[KD_VHL KD_CRBN];



% model fitting (at 1uM) and prediction (at 100nM) of Foretinib CRBN and VHL PROTACs

degradation=importdata('Foretinib_data.csv');  
degradation_data=degradation.data(:,1:4);
affinity_data=degradation.data(:,5:6)/1000;          % target KD; uM
halflife=degradation.data(:,7);               % target protein half-life; hr
copynumber=degradation.data(:,8);             % target protein copy number; 
kdp=log(2)/halflife/60;                    % target endogeneous protein degradation rate; min^(-1)
concentration=copynumber/V/NA_MOLE*10^6;   % target protein concentration; uM

dose_1=1;     % PROTAC dose; uM
dose_2=0.1;   % PROTAC dose; uM

T=24*60;   % timespan; min

options=odeset('events',@events_time);
tic;

kpr_series=10.^[-2:0.01:2];

min_kpr=zeros(length(halflife),2);
min_fc=zeros(length(halflife),4);

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

csvwrite('Foretinib_min_kpr.csv',min_kpr);
csvwrite('Foretinib_min_fc.csv',min_fc);



% model fitting (at 1uM) and prediction (at 100nM) of TAE684 CRBN PROTACs

TAE684=importdata('TAE684_data.csv');
degradation_data=TAE684.data(:,2);

affinity_data=TAE684.data(:,3)/1000;       % target KD; uM

halflife=TAE684.data(:,4);                 % target protein half-life; hr
copynumber=TAE684.data(:,5);               % target protein copy number; 
kdp=log(2)/halflife/60;                    % target endogeneous protein degradation rate; min^(-1)
concentration=copynumber/V/NA_MOLE*10^6;   % target protein concentration; uM

dose_1=1;     % PROTAC dose; uM
dose_2=0.1;   % PROTAC dose; uM

T1=5*60;      % timespan; min
T2=4*60;      % timespan; min

min_kpr=zeros(length(halflife),1);
min_fc=zeros(length(halflife),2);

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

csvwrite('TAE684_min_kpr.csv',min_kpr);
csvwrite('TAE684_min_fc.csv',min_fc);

