function [CssExtra]=ASET_ExtrapCss(V,S,Css,Cut_back,St,CssExtraM,P,C_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapole the different method possible in ASET. The methods
% are: None, constant value, Linear (3pt data), Rouse distribution.

% by Dominguez Ruben, L. FICH- UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surfsize=size(C_back.zsurfback,1);
bottomsize=size(C_back.zbottomback,1);

for j=1:St.mback
    if Cut_back.Backcut(1,j)<=6 %when you have least 6 valid data for ensemble
        if Cut_back.Backcut1(1,j)==0 
            CssExtra.surf(1:surfsize,j)=nan;
        else   
            for e=1:length(C_back.zsurfback(:,j));
                if isnan(C_back.zsurfback(e,j))
                    CssExtra.surf(e,j)=nan;
                else
                    CssExtra.surf(e,j)=Css(Cut_back.Backcut1(1,j),j); 
                end
            end                     
        end

        if Cut_back.Backcut1(1,j)==0
            CssExtra.bottom(1:bottomsize,j)=nan;
        else               

            for e=1:length(C_back.zbottomback(:,j));
                if isnan(C_back.zbottomback(e,j))
                    CssExtra.bottom(e,j)=nan;
                else
                    CssExtra.bottom(e,j)=Css(Cut_back.Backcut2(1,j),j);
                end
            end
        end

    else %when you have more than 6 valid data for ensemble
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sediment Extrapolation Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (strcmp(CssExtraM,'CON'))%Constany Extrapolation method
            [CssEx]=ASET_ExtrapCssConstant(Css,C_back,j,Cut_back);

            CssExtra.bottom(:,j)=CssEx.bottom; 

            CssExtra.surf(:,j)=CssEx.surf;

        elseif (strcmp(CssExtraM,'LSI'))

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %3pt Sediment Extrapolation Method 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [CssEx]=ASET_ExtrapCss3pt(V,Css,C_back,j,Cut_back);

            CssExtra.bottom(:,j)=CssEx.bottom; 
            CssExtra.predbottom(:,j)=CssEx.predbottom;

            CssExtra.surf(:,j)=CssEx.surf;
            CssExtra.predsurf(:,j)=CssEx.predsurf;

        elseif (strcmp(CssExtraM,'RSI'))

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Rouse Sediment Extrapolation Method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [CssEx]=ASET_ExtrapCssRouse(V,Css,C_back,P,S,j,Cut_back);

             CssExtra.pred(:,j)= CssEx.pred;
             CssExtra.bottom(:,j) = CssEx.bottom;
             CssExtra.surf(:,j) = CssEx.surf;
             CssExtra.Ca(:,j)=CssEx.Ca;
             CssExtra.NRouse(:,j)=CssEx.NRouse;
             CssExtra.a(:,j)=CssEx.a;
             CssExtra.zpred(:,j)=CssEx.zpred;
             CssExtra.acoezpred(:,j)=CssEx.acoezpred;
             CssExtra.acoezbottom(:,j)= CssEx.acoezbottom;
             CssExtra.acoezsurf(:,j) = CssEx.acoezsurf;
             CssExtra.rcuad(:,j)=CssEx.rcuad;
             CssExtra.ustarrouse(:,j)=CssEx.ustarrouse;

        end
        
        % If extrapolation is bad use the constant extrapolation method
        if isreal(CssExtra.bottom(:,j)) |  isreal(CssExtra.surf(:,j))
        else
            [CssEx]=ASET_ExtrapCssConstant(Css,C_back,j,Cut_back);
            CssExtra.bottom(:,j)=CssEx.bottom; 
            CssExtra.surf(:,j)=CssEx.surf;
        end
    end

            
end

%Control the value calculate and delete the bad data
CssExtra.bottom(CssExtra.bottom==Inf)=nan;
CssExtra.bottom(CssExtra.bottom==-Inf)=nan;
CssExtra.bottom(CssExtra.bottom<0)=nan;
CssExtra.surf(CssExtra.surf==Inf)=nan;
CssExtra.surf(CssExtra.surf==-Inf)=nan;
CssExtra.surf(CssExtra.surf<0)=nan;