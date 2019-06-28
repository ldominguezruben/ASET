function [C_vel]=ASET_ArrayResizeVelExtrap(V,St,Cut_vel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determinate the cell on the extrapolation zone. Its posible
% calculate in secton with variable cell size. For Velocity variable

%by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:St.mback
    for i=1:St.nback
        C_vel.zfrombottom(i,j)=V.mcsBed(1,j)-V.mcsDepth(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extrapolation zones
%Bottom
for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
         ncellsbottom(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut2(1,j),j)+...
             (St.cell(Cut_vel.Velcut2(1,j),j)/2)))/St.cell(Cut_vel.Velcut2(1,j),j);
    end
end

C_vel.zbottomvel=nan(max(round(ncellsbottom)+1),St.mvel);

for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
        
    % Bottom Location

        if V.mcsBed(1,j)-V.mcsDepth(Cut_vel.Velcut2(1,j),j)<St.cell(Cut_vel.Velcut2(1,j),j)
            C_vel.zbottomvel(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut2(1,j),j)))/2;%from bottom
        else
             ncellsbottom(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut2(1,j),j)+...
                (St.cell(Cut_vel.Velcut2(1,j),j)/2)))/St.cell(Cut_vel.Velcut2(1,j),j);
            if ncellsbottom(1,j)<0
                C_vel.zbottomvel(:,j)=nan;
            else
            for uu=1:fix(ncellsbottom(1,j))            
                C_vel.zbottomvel(uu,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut2(1,j),j)))-(uu*St.cell(Cut_vel.Velcut2(1,j),j));
            end

            if ncellsbottom(1,j)-floor(ncellsbottom(1,j))~=0
               C_vel.zbottomvel(uu+1,j)=C_vel.zbottomvel(uu,j);%
                if C_vel.zbottomvel(uu+1,j)<V.mcsBed(1,j)%
                   C_vel.zbottomvel(uu+1,j)=nan;%
                else
                end
            else
            end
            C_vel.zbottomvel(uu+2,j)=0;%last data on 5%H
            end
        end
    end
end

C_vel.zbottomvel(C_vel.zbottomvel==0)=nan;%nunca puede ser 0

%Prediction vector bottom 
for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
        C_vel.zbottom(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut2(1,j),j)))/2;%from bottom
        C_vel.zpredbottomvel(:,j)=linspace(C_vel.zbottom(1,j)-1,C_vel.zbottom(1,j)+1,20);
    end
end

%Determinate the cellsize of bottom unmeasured zone
for j=1:St.mvel
    for uu=1:size(C_vel.zbottomvel,1)  
        if uu==1
            C_vel.cellbottomVel(uu,j)=abs(C_vel.zbottomvel(uu+1,j)-C_vel.zbottomvel(uu,j));
        else
            C_vel.cellbottomVel(uu,j)=abs(C_vel.zbottomvel(uu,j)-C_vel.zbottomvel(uu-1,j));%saber que la mitad de la celda no esta considerada
        end
    end
end

%%
%Position surface points

for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
        ncellssurf(1,j)=(V.mcsDepth(Cut_vel.Velcut1(1,j),j)-(St.cell(1,j)/2))/St.cell(1,j);
    end
end

C_vel.zsurfvel=nan(max(round(ncellssurf)),St.mvel);

for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
        if V.mcsDepth(Cut_vel.Velcut1(1,j),j)<St.cell(1,j)%only for RiverRay ADCP
            C_vel.zsurfvel(1,j)=V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut1(1,j),j)/2);
        else% for the rest of ADCP
            ncellssurf(1,j)=(V.mcsDepth(Cut_vel.Velcut1(1,j),j)-(St.cell(1,j)/2))/St.cell(1,j);
            for uu=1:fix(ncellssurf(1,j))            
                C_vel.zsurfvel(uu,j)=uu*St.cell(1,j)+(V.mcsBed(1,j)-V.mcsDepth(Cut_vel.Velcut1(1,j),j));
            end
            if ncellssurf(1,j)-floor(ncellssurf(1,j))~=0
               C_vel.zsurfvel(uu+1,j)=(V.mcsBed(1,j)-C_vel.zsurfvel(uu,j))/2+C_vel.zsurfvel(uu,j);
            else
                
            end
        end
    else
    end
    C_vel.zsurfvel(:,j)=flipud(C_vel.zsurfvel(:,j));
end

%Delete bad data
C_vel.zsurfvel(C_vel.zsurfvel==0)=nan;

%prediction vector bottom points
for j=1:St.mvel
    if Cut_vel.Velcut(1,j)>0
        C_vel.zsurf(1,j)=V.mcsBed(1,j)-(V.mcsDepth(Cut_vel.Velcut1(1,j),j)/2);
        C_vel.zpredsurfvel(:,j)=linspace(C_vel.zsurf(1,j)-1,C_vel.zsurf(1,j)+1,20);
    end
end

%Determinate the cellsize of bottom unmeasured zone

C_vel.cellsurfVel=abs(diff(C_vel.zsurfvel,1));
t=size(C_vel.cellsurfVel,1);
if t<2
else
    C_vel.cellsurfVel(t,:)=nanmax(C_vel.cellsurfVel(1:t-1,:));%
    C_vel.cellsurfVel(t+1,:)=nanmax(C_vel.cellsurfVel);%add one for the last value
end
end