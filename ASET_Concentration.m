function [Css,Scb,Cal,P]=ASET_Concentration(M,V,S,Cut_back,R,St,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculate the concetration using the Deines 1999 equation,
%and Szupiany et al., 2018 correction. Also do the calibration of signal
%with concetration.

%by Dominguez Ruben, UNL, Argentina.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start code
if M==1 %Not Calibration. Transport and Monte Carlo
%%
%Correlation equation (Deines, 1999)
%SCB=10log10(Ms)=
%kc*E-Cr+2*(alphaw+alphas)*R+10log10((Tt*R^2*psi^2)/(L*Pt))+C-10log10(ks^2)
%        First Term          Second Term                      third term 
                          
T_2=2*P.Ac.ALPHAw*P.Ac.r;%attenuation due water iscosity
T_3=10*log10((273+P.Ph.Temp)*(P.Ac.r.^2).*(P.Ac.psi.^2)./(P.El.Pt*P.El.L));%Second term
T_4=10*log10(P.Ph.ks);%third term

%iteration depends of frequencies
if P.Ac.FREQ==1228.8;
    for j=1:St.mback
        for i=1:St.nback
            if isnan(V.mcsBack(i,j)) | V.mcsBack(i,j)==255
            Scb(i,j)=nan;
            P.Ac.ALPHAg(i,j)=nan;
            Css(i,j)=nan;     
            
            elseif i==Cut_back.Backcut1(j)
                
                for e=1:50%near field correction see Dominguez Ruben ASET paper first cell     
                   if e==1;                    
                    Scb(i,j)=(V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAg(i,j))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;   
                    Css(i,j)=10^(Scb(i,j)/10);% kg/m3

                    % Fitting  coarse sediment
                    ALPHAg1=(10*log10(exp(1)^2));
                    ALPHAg2=((P.Ac.k*(Css(i,j)*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                    ALPHAg3=((0.4*(Css(i,j)*1000))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                    P.Ac.ALPHAgc(e)=ALPHAg1*(ALPHAg2+ALPHAg3);

                   elseif e>1;

                    Scb(i,j)=(V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAgc(e-1))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;  
                    Css(i,j)=10^(Scb(i,j)/10);%kg/m3

                    % Fitting  coarse sediment
                    ALPHAg1=(10*log10(exp(1)^2));
                    ALPHAg2=((P.Ac.k*(Css(i,j)*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                    ALPHAg3=((0.4*(Css(i,j)*1000))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                    P.Ac.ALPHAgc(e)=ALPHAg1*(ALPHAg2+ALPHAg3);

                   end       

                end
            
            
            P.Ac.ALPHAg(i,j)=P.Ac.ALPHAgc(e);
            
            elseif i>Cut_back.Backcut1(j)
                             
                Scb(i,j)=((V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr)+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAg(i-1,j))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;
                Css(i,j)=10^(Scb(i,j)/10);%kg/m3

                % Fixing  coarse sediment
                ALPHAg1=(10*log10(exp(1)^2));
                ALPHAg2=((P.Ac.k*(nanmean(Css(Cut_back.Backcut1(j):i,j))*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                ALPHAg3=((0.4*(nanmean(Css(Cut_back.Backcut1(j):i,j))*1000))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                P.Ac.ALPHAg(i,j)=ALPHAg1*(ALPHAg2+ALPHAg3); 
                                
            end
        end
    end
    
            elseif P.Ac.FREQ==614;
                                
            for j=1:St.mback
            for i=1:St.nback
            if V.mcsBack(i,j)==255;
                Scb(i,j)=nan;
                P.Ac.ALPHAg(i,j)=nan;
                Css(i,j)=nan;
            elseif i==Cut_back.Backcut1(j)
                
                Scb(i,j)=((V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr)+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAg(i,j))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;
                Css(i,j)=10^(Scb(i,j)/10);%kg/m3  

            elseif i==Cut_back.Backcut1(j)+1
            for e=1:50%correction of alhphas2 
               if e==1;
                    
                Scb(i,j)=(V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAg(i,j))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;  
                Css(i,j)=10^(Scb(i,j)/10);% kg/m3

                % Fixing  coarse sediment
                ALPHAg1=(10*log10(exp(1)^2));
                ALPHAg2=((P.Ac.k*(Css(i,j)*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                ALPHAg3=((0.4*(Css(i,j)*1000))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                P.Ac.ALPHAgc(e)=ALPHAg1*(ALPHAg2+ALPHAg3);
            
               elseif e>1;
                   
                Scb(i,j)=(V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAgc(e-1))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;  
                Css(i,j)=10^(Scb(i,j)/10);%kg/m3

                % Fixing  coarse sediment
                ALPHAg1=(10*log10(exp(1)^2));
                ALPHAg2=((P.Ac.k*(Css(i,j)*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                ALPHAg3=((0.4*(Css(i,j)*1000))/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                P.Ac.ALPHAgc(e)=ALPHAg1*(ALPHAg2+ALPHAg3);
                   
               end    
            
            end           
            
            P.Ac.ALPHAg(i,j)=P.Ac.ALPHAgc(e);
            
            elseif i>Cut_back.Backcut1(j)+1
                
                Scb(i,j)=((V.mcsBack(i,j).*P.Ac.mean_beam)-P.Ac.kcEr)+T_2(i,j)+2*(P.Ac.ALPHAf+P.Ac.ALPHAg(i-1,j))*P.Ac.r(i,j)+T_3(i,j)+P.Ac.Ccoef-T_4;
                Css(i,j)=10^(Scb(i,j)/10);%kg/m3

                % Fixing  coarse sediment 
                ALPHAg1=(10*log10(exp(1)^2));
                ALPHAg2=((P.Ac.k*(nanmean(Css(Cut_back.Backcut1(j)+1:i,j))*1000)*((P.Ph.s-1)^2))/(2*P.Ph.Rhos))*(P.Ac.sc/((P.Ac.sc^2)+((P.Ph.s+P.Ac.tauc)^2)));
                ALPHAg3=((0.4*(nanmean(Css(Cut_back.Backcut1(j)+1:i,j)))*1000)/(S.Dc*P.Ph.Rhos))*(((P.Ac.k*S.Dc/2)^4)/(1+(1.3*(P.Ac.k*S.Dc/2)^2)+(0.24*(P.Ac.k*S.Dc/2)^4)));

                P.Ac.ALPHAg(i,j)=ALPHAg1*(ALPHAg2+ALPHAg3);
                
                                 
            end
            end
        end
end

kcEr=[];%calculate not calibration
Cal=[];%calibrate not calibration

elseif M==2%Calibration 
    
    %%
    %Correlation equation (Deines, 1999)
    %SCB=10log10(Ms)=kc*E-Cr+2*(alphaw+alphas)*R+10log10((Tt*R^2*psi^2)/(L*Pt))+C-10log10(ks^2)
    for i=1:size(R,2)

        %%
        %Intensity array
        ave=[R{i}.Back.Beam1(:,:,1).*P{i}.Ac.beam_1 R{i}.Back.Beam2(:,:,1).*P{i}.Ac.beam_2 R{i}.Back.Beam3(:,:,1).*P{i}.Ac.beam_3 R{i}.Back.Beam4(:,:,1).*P{i}.Ac.beam_4];
        %Transforma los cero en nan
        ave(ave==0)=nan;
        R{i}.Back.BeamAve=nanmean(ave,2);       
        
        %Delete temporaly vars
        clear ave
        
        %Correlation equation (Deines, 1999)
        %SCB=10log10(Ms)=
        %kc*E-Cr+2*(alphaw+alphas)*R+10log10((Tt*R^2*psi^2)/(L*Pt))+C-10log10(ks^2)
        %        First Term          Second Term                      third term 
    
        T_2{i}=2*P{i}.Ac.ALPHAw.*P{i}.Ac.r;%First term
        T_3{i}=10*log10((273+P{i}.Ph.Temp)*(P{i}.Ac.r.^2).*(P{i}.Ac.psi.^2)./(P{i}.El.Pt*P{i}.El.L));%Second term
        T_4{i}=10*log10(P{i}.Ph.ks);%third term

    %iteration depends of frequencies
        for ii=1:St{i}.nback
            if isnan(R{i}.Back.BeamAve(ii,1))
                P{i}.Ac.ALPHAg(ii,1)=nan;  
                Cal.kcEr{i}(ii,1)=nan;
                Cal.FCB{i}(ii,1)=nan;
                Cal.kE{i}(ii,1)=nan;
                Cal.ST{i}(ii,1)=nan;
            else                
                Cal.kcEr{i}(ii,1)=R{i}.Back.BeamAve(ii,1)-(10*log10(S.Css(i)/1000)-T_2{i}(ii,1)-2*(P{i}.Ac.ALPHAf+P{i}.Ac.ALPHAg(ii,1))*P{i}.Ac.r(ii,1)-...
                    T_3{i}(ii,1)-P{i}.Ac.Ccoef+T_4{i});                                         
                Cal.FCB{i}(ii,1)=20*log10(P{i}.Ac.psi(ii,1)*P{i}.Ac.r(ii,1))+T_2{i}(ii,1)+R{i}.Back.BeamAve(ii,1);        
                Cal.kE{i}(ii,1)=10*log10(S.Css(i)/1000)-T_2{i}(ii,1)-2*(P{i}.Ac.ALPHAf+P{i}.Ac.ALPHAg(ii,1))*P{i}.Ac.r(ii,1)-T_3{i}(ii,1)-P{i}.Ac.Ccoef+T_4{i};
                Cal.ST{i}(ii,1)=Cal.FCB{i}(ii,1)+2*(P{i}.Ac.ALPHAf+P{i}.Ac.ALPHAg(ii,1))*P{i}.Ac.r(ii,1);                                          
            end
                Cal.ST_mean{i}=nanmean(Cal.ST{i}(:,1));
        end

           [minimo,Cal.indice(i)]=min(abs(R{i}.Depth_m(:,1)-S.depthdata(i)));
            ave=[R{i}.Back.Beam1(:,:,1) R{i}.Back.Beam2(:,:,1) R{i}.Back.Beam3(:,:,1) R{i}.Back.Beam4(:,:,1)];
            m(i)=nansum(isfinite(ave(Cal.indice(i),:)));

           if (m(i)/size(ave,2))*100>10% Control minimum 10% of data is incoporate in the cells close to the bottom
                Cal.kcEr_puntual(i,1)=Cal.kcEr{i}(Cal.indice(i),1);
           end

        clear ave
    end

    Cal.kcEr_mean=nanmean(Cal.kcEr_puntual);

    %%
    for i=1:size(R,2)
        for ii=1:St{i}.nback
            if isnan(R{i}.Back.BeamAve(ii,1))
                P{i}.Ac.ALPHAg(ii,1)=nan;     
                Cal.FCB1{i}(ii,1)=nan;
                Cal.SCB{i}(ii,1)=nan;
                Cal.Css{i}(ii,1)=nan;
            else 
                Cal.FCB1{i}(ii,1)=T_3{i}(ii,1)+T_2{i}(ii,1)+R{i}.Back.BeamAve(ii,1);
                Cal.SCB{i}(ii,1)=Cal.FCB1{i}(ii,1)+2*(P{i}.Ac.ALPHAf+P{i}.Ac.ALPHAg(ii,1))*P{i}.Ac.r(ii,1)-T_4{i}+P{i}.Ac.Ccoef-Cal.kcEr_mean;
                Cal.Css{i}(ii,1)=10^(Cal.SCB{i}(ii,1)/10);          
            end
        end
    end


Css=[];%not calculate calibration
Scb=[];%not calculate calibration
end



    
    