function [VelExtra]=ASET_ExtrapV3pt(V,C_vel,j,Cut_vel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapolate in surface and bottom the velocity value using
% the 3 first and the last valid data.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linear Velocity Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
%Surface Extrapolation

Velsurf=V.mcsMag(Cut_vel.Velcut1(1,j):Cut_vel.Velcut1(1,j)+2,j)./100;%m/s
Depthsurf=V.mcsBed(1,j)-V.mcsDepth(Cut_vel.Velcut1(1,j):Cut_vel.Velcut1(1,j)+2,j);

linearsurf=polyfit(Depthsurf,Velsurf,1);

%%%%%

%Define vector predection          
VelExtra.predsurf=polyval(linearsurf,C_vel.zpredsurfvel(:,j));

VelExtra.surf=polyval(linearsurf,C_vel.zsurfvel(:,j));

%%Bottom extrapolation

Velbottom=V.mcsMag(Cut_vel.Velcut2(1,j)-2:Cut_vel.Velcut2(1,j),j)./100;
Depthbottom=V.mcsBed(1,j)-V.mcsDepth(Cut_vel.Velcut2(1,j)-2:Cut_vel.Velcut2(1,j),j);

linearbottom=polyfit(Depthbottom,Velbottom,1);

%Define vector predection  
VelExtra.predbottom=polyval(linearbottom,C_vel.zpredbottomvel(:,j));

VelExtra.bottom=polyval(linearbottom,C_vel.zbottomvel(:,j));

    