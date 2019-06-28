function [VelExtra]=ASET_ExtrapVLawOfTheWall(V,C_vel,j,Cut_vel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the extrapolate data using the Law of the Wall
% close to the surface and bottom.
% LPVI

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Law of the Wall Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHOD LPVI  vel=a*ln(z)+b Equation

%Input variables 
depthlog=log(C_vel.zfrombottom(Cut_vel.Velcut1(1,j):Cut_vel.Velcut2(1,j),j));%m
u=V.mcsMag(Cut_vel.Velcut1(1,j):Cut_vel.Velcut2(1,j),j)./100;%m/s
VelExtra.zpred = linspace(0,V.mcsBed(1,j),100);%m

 
%Extrapolation method

Velfit=polyfit(depthlog,u,1);
Vely=Velfit(1)*depthlog+Velfit(2);

VelExtra.pred = polyval(Velfit,log(VelExtra.zpred));%For graphics
VelExtra.bottom = polyval(Velfit,log(C_vel.zbottomvel(:,j)));%Bottom extrapolation
VelExtra.surf = polyval(Velfit,log(C_vel.zsurfvel(:,j)));%Surface extrapolation

VelExtra.rsqvel = max(0,1 - sum((u(:)-Vely(:)).^2)/sum((u(:)-mean(u(:))).^2));


kappa=0.41;%constant von Karman

%Variables
VelExtra.ustarwl=Velfit(1)*kappa;%m/s
VelExtra.ordenada=Velfit(2);

clear nanIndex
