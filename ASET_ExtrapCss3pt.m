function [CssExtra]=ASET_ExtrapCss3pt(V,Css,C_back,j,Cut_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapolate the linear methods. This methods use the
% similar relation of velocity. Ask about this methos

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linear Sediment Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%surface sediment extrapolation

Csssurf=Css(Cut_back.Backcut1(1,j):Cut_back.Backcut1(1,j)+2,j);
Depthsurf=V.mcsBed(1,j)-V.mcsDepth(Cut_back.Backcut1(1,j):Cut_back.Backcut1(1,j)+2,j);

linearsurf=polyfit(Depthsurf,Csssurf,1);

CssExtra.surf=polyval(linearsurf,C_back.zsurfback(:,j));

%Prediction
CssExtra.predsurf=polyval(linearsurf,C_back.zpredsurfback(:,j));

for e=1:length(C_back.zsurfback(:,j));
    if CssExtra.surf(e,1)<0
        CssExtra.surf(e,1)=Css(Cut_back.Backcut1(1,j),j)/100;
    else   
    end
end
    

%Bottom sediment extrapolation

Cssbottom=Css(Cut_back.Backcut2(1,j)-2:Cut_back.Backcut2(1,j),j);
Depthbottom=V.mcsBed(1,j)-V.mcsDepth(Cut_back.Backcut2(1,j)-2:Cut_back.Backcut2(1,j),j);

linearbottom=polyfit(Depthbottom,Cssbottom,1);     

CssExtra.bottom=polyval(linearbottom,C_back.zbottomback(:,j)); 

%Prediction
CssExtra.predbottom=polyval(linearbottom,C_back.zpredbottomback(:,j));