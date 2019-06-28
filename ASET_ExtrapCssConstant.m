function [CssExtra]=ASET_ExtrapCssConstant(Css,C_back,j,Cut_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapolate the Ms2 using the first and last date in the 
% surface and bottom zone.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Surface zone Ms2 concentration
for e=1:length(C_back.zsurfback(:,j));
    if isnan(C_back.zsurfback(e,j))
        CssExtra.surf(e,1)=nan;
    else
        CssExtra.surf(e,1)=Css(Cut_back.Backcut1(1,j),j);       
    end
end

%Bottom zone Ms2 concetration
for e=1:length(C_back.zbottomback(:,j));
    if isnan(C_back.zbottomback(e,j))
        CssExtra.bottom(e,1)=nan;
    else
        CssExtra.bottom(e,1)=Css(Cut_back.Backcut2(1,j),j);
    end
end   

