function [CssExtra]=ASET_ExtrapCssRouse(V,Css,C_back,P,S,j,Cut_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapolate the coarse Concetration using the Rouse Distribution.
% See Rouse (1937). RoD

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surfsize=size(C_back.zsurfback,1);
bottomsize=size(C_back.zbottomback,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rouse Sediment Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CssExtra.C=Css(Cut_back.Backcut1(1,j):Cut_back.Backcut2(1,j),j);%kg/m3 from bottom

%Indexes
kappa=0.41;%von Karman constant
H=V.mcsBed(1,j);%Depth per each ensemble
CssExtra.a=0.05*H;%Ms2 reference from 5% from Bottomn for Rouse Distribution (Rouse, 1937)

%z values (input variables)
CssExtra.zpred = linspace(0.02*V.mcsBed(1,j),0.98*V.mcsBed(1,j),100); 
CssExtra.acoe=((H-C_back.zfrombottom(1:Cut_back.Backcut2(1,j),j))./C_back.zfrombottom(1:Cut_back.Backcut2(1,j),j))./((H-CssExtra.a)/CssExtra.a);                
CssExtra.acoezpred=(((H-CssExtra.zpred)./CssExtra.zpred)./((H-CssExtra.a)/CssExtra.a));
CssExtra.acoezbottom=(((H-C_back.zbottomback(:,j))./C_back.zbottomback(:,j))./((H-CssExtra.a)/CssExtra.a));
CssExtra.acoezsurf=(((H-C_back.zsurfback(:,j))./C_back.zsurfback(:,j))./((H-CssExtra.a)/CssExtra.a));

%% Fit of concentration 
if P.Ac.FREQ==1228.8;

    [Sedfit]=polyfit(log10(CssExtra.acoe),log10(CssExtra.C),1);
    Sedy=Sedfit(1)*log10(CssExtra.acoe)+Sedfit(2);

    CssExtra.pred =  polyval(Sedfit,log10(CssExtra.acoezpred));
    CssExtra.bottom = polyval(Sedfit,log10(CssExtra.acoezbottom));
    CssExtra.surf = polyval(Sedfit,log10(CssExtra.acoezsurf));

    CssExtra.rcuad = max(0,1 - sum((log10(CssExtra.C(:))-Sedy(:)).^2)/sum((log10(CssExtra.C(:))-mean(log10(CssExtra.C(:)))).^2));

    %Conversion         
    CssExtra.bottom(1:bottomsize,1)=10.^(CssExtra.bottom);
    CssExtra.surf(1:surfsize,1)=10.^(CssExtra.surf);
    CssExtra.pred=10.^(CssExtra.pred);
    CssExtra.Ca = 10^Sedfit(2);
    CssExtra.NRouse=10^Sedfit(1);

    %%%%%%%%

    if Sedfit(1)>0
    elseif Sedfit(1)<0
        for t=1:surfsize
            CssExtra.surf(t,1)=nanmean(Css(Cut_back.Backcut1(1,j):Cut_back.Backcut1(1,j)+2,j));
        end
        for e=1:bottomsize
            if isnan(CssExtra.acoezbottom(e,1))
            else
                CssExtra.bottom(e,1)=nanmean(Css(Cut_back.Backcut2(1,j)-2:Cut_back.Backcut2(1,j),j));
                CssExtra.Ca = nan;
                CssExtra.NRouse= nan;
            end
        end
    end

elseif P.Ac.FREQ==614;

     Sedfit=polyfit(log10(CssExtra.acoe(2:end,1)),log10(CssExtra.C(2:end,1)),1);
     Sedy=Sedfit(1)*log10(CssExtra.acoe(2:end,1))+Sedfit(2);

     CssExtra.pred= polyval(Sedfit,log10(CssExtra.acoezpred));
     CssExtra.bottom = polyval(Sedfit,log10(CssExtra.acoezbottom));
     CssExtra.surf=polyval(Sedfit,log10(CssExtra.acoezsurf));

     CssExtra.rcuad = max(0,1 - sum((log10(CssExtra.C(2:end,1))-Sedy(:)).^2)/sum((log10(CssExtra.C(2:end,1))-mean(log10(CssExtra.C(2:end,1)))).^2));

    CssExtra.bottom(1:bottomsize,1)=10.^(CssExtra.bottom);
    CssExtra.surf(1:surfsize,1)=10.^(CssExtra.surf);
    CssExtra.pred=10.^(CssExtra.pred);
    CssExtra.Ca = 10^Sedfit(2);
    CssExtra.NRouse=10^Sedfit(1);


% Case that the extrapolation downstream per ensemble
    if Sedfit(1)>0
    elseif Sedfit(1)<0
        for t=1:surfsize
            CssExtra.surf(t,1)=nanmean(Css(Cut_back.Backcut1(1,j):Cut_back.Backcut1(1,j)+2,j));
        end
        for e=1:bottomsize
        if isnan(CssExtra.acoezbottom(e,1))
            CssExtra.bottom(e,1)=nan;
        else
            CssExtra.bottom(e,1)=nanmean(Css(Cut_back.Backcut2(1,j)-2:Cut_back.Backcut2(1,j),j));
            CssExtra.Ca = nan;
            CssExtra.NRouse= nan;
        end
        end

    end
end
       
%%
% Sediment Parameters 
g=9.8;%m2/s              

%shear velocity calculate
kappa=0.41;%adimensional

%Settling velocity equations
if S.Dc<10E-4
    CssExtra.wsettling=((P.Ph.s-1)*g*S.Dc^2)/(18*P.Ph.nu);
elseif 10E-4<=S.Dc<=10E-3
    CssExtra.wsettling=(10*P.Ph.nu)/(S.Dc)*((1+0.01*(P.Ph.s-1)*g*S.Dc^3)/(18*P.Ph.nu^2)-1)^0.5;
elseif S.Dc>10E-3
    CssExtra.wsettling=1.1*((P.Ph.s-1)*g*S.Dc)^0.5;
end 
%ustar by Rouse
CssExtra.ustarrouse=CssExtra.wsettling/(kappa.*CssExtra.NRouse);%m s-1