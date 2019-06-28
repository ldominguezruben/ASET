function [P]=ASET_Parameters(M,V,S,R,Inst,Cfg,Sensor,St,ADCPdata,Mont,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This funtion calculate the basic parameters to calculate: Acoustic,
%physics, and electric.

%Dominguez Ruben, L. UNL (03/2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M==1%single file or Multiselecte one by one
    Depth_m=V.mcsDepth;
elseif M==2%calibration
    Depth_m=R{i}.Depth_m;
end
    
%Temperature parameters [c]
Sensor.temperature_degc(Sensor.temperature_degc==-32768)=[];%delete empty cell
P.Ph.Temp=nanmean(Sensor.temperature_degc);

%Parameter Acoustic and Physics
%Frequency ADCP
P.Ac.FREQ=Inst.freq(1);

if P.Ac.FREQ==600
    P.Ac.FREQ=614;%Frenquency equipment
    P.Ac.Draft_ADCP=ADCPdata.R.draft(i);
    P.Ac.ALPHAt=0.0365;%ALPHA coefficient
    P.Ac.beam_1=ADCPdata.rssi_beam1(i);
    P.Ac.beam_2=ADCPdata.rssi_beam2(i);
    P.Ac.beam_3=ADCPdata.rssi_beam3(i);
    P.Ac.beam_4=ADCPdata.rssi_beam4(i);
    beam=[P.Ac.beam_1;P.Ac.beam_2;P.Ac.beam_3;P.Ac.beam_4];
    P.Ac.mean_beam=mean(beam);%average beam
    P.Ac.kcEr=ADCPdata.kcEr;
    P.Ac.Ccoef=str2double(ADCPdata.ccoef);
elseif P.Ac.FREQ==1200
    P.Ac.FREQ=1228.8;%frequency equipment
    P.Ac.Draft_ADCP=ADCPdata.R.draft(i);%boat distance [m]
    P.Ac.ALPHAt=0.0269875;%ALPHA coefficient
    P.Ac.beam_1=ADCPdata.rssi_beam1(i);
    P.Ac.beam_2=ADCPdata.rssi_beam2(i);
    P.Ac.beam_3=ADCPdata.rssi_beam3(i);
    P.Ac.beam_4=ADCPdata.rssi_beam4(i);
    beam=[P.Ac.beam_1;P.Ac.beam_2;P.Ac.beam_3;P.Ac.beam_4];
    P.Ac.mean_beam=mean(beam);%average beam
    P.Ac.kcEr=ADCPdata.kcEr;
    P.Ac.Ccoef=str2double(ADCPdata.ccoef);
end

%Sound parameters
P.Ac.vson=1500;%speed of sound [m s-1]
P.Ac.lambda=P.Ac.vson/(P.Ac.FREQ*1000);%wave longitud
P.Ac.k=2*pi/P.Ac.lambda;%k coefficient

%Density material and variables
P.Ph.Rhos=2650000;%sediment density (cuarzo)[mg l-1]
P.Ph.x=P.Ac.k*(S.Dc(i)/2);%dimension coefficient
P.Ph.f=1.25*(P.Ph.x.^2);%form factor
P.Ph.ks=(3/(16*pi*P.Ph.Rhos/1000))*((P.Ph.f.^2)/(S.Dc(i)/2));%ks^2

%%
  
P.Ph.Rhof=1000*((999.83952+16.945176*P.Ph.Temp-7.9870401E-3*P.Ph.Temp^2+46.170461E-6*P.Ph.Temp^3+...
    105.56302E-9*P.Ph.Temp^4-280.54253E-12*P.Ph.Temp^5)/(1+16.89785E-3*P.Ph.Temp));% (Kell, 1975) valid from 0ºC to 150ºC%. Water Density [mg/l]

P.Ph.nu=1.79E-6/(1+0.03368*P.Ph.Temp+0.00021*P.Ph.Temp^2);%Kinematic Viscosity water [m2/s] -Garcia (2008)

%%
%ALPHA attenuation variables 

%ALPHAw fluid 
P.Ac.ALPHAw=(10*log10(exp(1)^2)*(0.00000338*(P.Ac.FREQ^2)))/(21.9*(10^(6-(1520/(P.Ph.Temp+273))))); %see Latosisnski et al. (2014)

%ALPHAs fine and coarse material
P.Ph.s=P.Ph.Rhos/P.Ph.Rhof;
P.Ac.gama=(pi*P.Ac.FREQ*1000/P.Ph.nu)^0.5; %frequency seconds

P.Ac.sf=(9/(2*P.Ac.gama*S.Df(i)))*(1+(2/(P.Ac.gama*S.Df(i))));%add units
P.Ac.sc=(9/(2*P.Ac.gama*S.Dc(i)))*(1+(2/(P.Ac.gama*S.Dc(i))));%add units
P.Ac.tauf=0.5+(9./(2*P.Ac.gama*S.Df(i)));%tau fine material
P.Ac.tauc=0.5+(9./(2*P.Ac.gama*S.Dc(i)));%tau coarse material 

if M==2 %for calibration only Szupiany et al. (2019)
       
    %ALPHA fine material
    P.Ac.ALPHAf1=(10*log10(exp(1)^2));
    P.Ac.ALPHAf2=((P.Ac.k*(nanmean(S.Csf(1:i)))*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sf/((P.Ac.sf.^2)+((P.Ph.s+P.Ac.tauf).^2)));
    P.Ac.ALPHAf3=((0.4*(nanmean(S.Csf(1:i))))/(S.Df(i).*P.Ph.Rhos))*(((P.Ac.k*S.Df(i)/2).^4)./(1+(1.3*(P.Ac.k*S.Df(i)/2).^2)+(0.24*(P.Ac.k*S.Df(i)/2).^4)));
    P.Ac.ALPHAf=P.Ac.ALPHAf1*(P.Ac.ALPHAf2+P.Ac.ALPHAf3);%Alpha_s1 attenuation parameter from viscosity and scattering effect (fine material)

    %ALPHA coarse material
    P.Ac.ALPHAg1=(10*log10(exp(1)^2));
    P.Ac.ALPHAg2=((P.Ac.k*(nanmean(S.Css(1:i)))*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc.^2)+((P.Ph.s+P.Ac.tauc).^2)));
    P.Ac.ALPHAg3=((0.4*(nanmean(S.Css(1:i))))/(S.Dc(i).*P.Ph.Rhos))*(((P.Ac.k*S.Dc(i)/2).^4)./(1+(1.3*(P.Ac.k*S.Dc(i)/2).^2)+(0.24*(P.Ac.k*S.Dc(i)/2).^4)));
    P.Ac.ALPHAgr=P.Ac.ALPHAg1*(P.Ac.ALPHAg2+P.Ac.ALPHAg3);%Alpha_s2 attenuation parameter from viscosity and scattering effect (coarse material)

elseif M==1 %for the others modules

    %ALPHA fine material
    P.Ac.ALPHAf1=(10*log10(exp(1)^2));
    P.Ac.ALPHAf2=((P.Ac.k*((S.Csf))*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sf/((P.Ac.sf.^2)+((P.Ph.s+P.Ac.tauf).^2)));
    P.Ac.ALPHAf3=((0.4*((S.Csf)))/(S.Df.*P.Ph.Rhos))*(((P.Ac.k*S.Df/2).^4)./(1+(1.3*(P.Ac.k*S.Df/2).^2)+(0.24*(P.Ac.k*S.Df/2).^4)));
    P.Ac.ALPHAf=P.Ac.ALPHAf1*(P.Ac.ALPHAf2+P.Ac.ALPHAf3);%Alpha_s1 attenuation parameter from viscosity and scattering effect (fine material)

    %ALPHA coarse material
    P.Ac.ALPHAg1=(10*log10(exp(1)^2));
    P.Ac.ALPHAg2=((P.Ac.k*((S.Css))*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
    P.Ac.ALPHAg3=((0.4*((S.Css)))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));
    P.Ac.ALPHAgr=P.Ac.ALPHAg1*(P.Ac.ALPHAg2+P.Ac.ALPHAg3);%Alpha_s2 attenuation parameter from viscosity and scattering effect (coarse material)
end
 
P.Ac.ALPHAg = ones(St.nback,St.mback).*P.Ac.ALPHAgr;%array coarse ALPHA_S attenuation parameter

% Variables define SCB value
P.Ac.r = (Depth_m-P.Ac.Draft_ADCP)./cos((pi*20)/180);%projection 20 degrees beam depth [m]
P.Ac.R = P.Ac.r+((St.cell/4)/cos((pi*20)/180));%[m] correction phi
P.Ac.z = P.Ac.R*P.Ac.lambda/(pi*(P.Ac.ALPHAt^2));%

if M==2
    P.Ac.psi = (1+(1.35*P.Ac.z)+((2.5*P.Ac.z).^3.2))./((1.35*P.Ac.z)+((2.5*P.Ac.z).^3.2));%
else
    P.Ac.psi = ones(St.nback,St.mback);% 
end
    
%%
%Electric Parameter 
%Read .PDO variables file 
P.El.TC=nanmean(Sensor.xmitCurrent);%Transmit current

if P.El.TC==0 | isempty(P.El.TC)%Empty transmit current data. Put default value
    P.El.TC=50.5;
end

P.El.TV=nanmean(Sensor.xmitVoltage);%Transmit voltage

P.El.L=nanmean(Cfg.xmitPulse_cm)/100;%Transmit pulse
%variables depend of FREQ
if P.Ac.FREQ==1228.8;
    A1=0.011451;
    B1=0.253765;
elseif P.Ac.FREQ==614;
    A1=0.011451;
    B1=0.380667;
end

P.El.Pt=(P.El.TC*A1)*(P.El.TV*B1);%add definition

clear St{i}.nback St{i}.mback Depth_m
