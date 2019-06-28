function [C_back]=ASET_ArrayResizeBackExtrap(V,St,Cut_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determinate the cell on the extrapolation zone. Its posible
% calculate in secton with variable cell size. For Intensity

%by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:St.mvel
    for i=1:St.nvel
        C_back.zfrombottom(i,j)=V.mcsBed(1,j)-V.mcsDepth(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Build the bottom region

for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
        ncellsbottom(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut2(1,j),j)+(St.cell(Cut_back.Backcut2(1,j),j)/2)))/St.cell(Cut_back.Backcut2(1,j),j);
    end
end

C_back.zbottomback=nan(max(round(ncellsbottom)+1),St.mback);

for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
%Position bottom
        if V.mcsBed(1,j)-V.mcsDepth(Cut_back.Backcut2(1,j),j)<St.cell(Cut_back.Backcut2(1,j),j)
            C_back.zbottomback(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut2(1,j),j)))/2;%from bottom
        else
            ncellsbottom(1,j)=(0.95*V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut2(1,j),j)+(St.cell(Cut_back.Backcut2(1,j),j)/2)))/St.cell(Cut_back.Backcut2(1,j),j);

            if ncellsbottom(1,j)<0
                 C_back.zbottomback(:,j)=nan;
            else
            for uu=1:fix(ncellsbottom(1,j))            
                C_back.zbottomback(uu,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut2(1,j),j)))-uu*St.cell(Cut_back.Backcut2(1,j),j);
            end
            if ncellsbottom(1,j)-floor(ncellsbottom(1,j))~=0
                if Cut_back.Backcut2(1,j)+fix(ncellsbottom(1,j))<size(V.mcsDepth,1)
                     C_back.zbottomback(uu+1,j)=(V.mcsBed(1,j)-V.mcsDepth(Cut_back.Backcut2(1,j)+fix(ncellsbottom(1,j)),j))-((ncellsbottom(1,j)-floor(ncellsbottom(1,j)))/2*St.cell(end,j));
                else
                    C_back.zbottomback(uu+1,j)=(V.mcsBed(1,j)-V.mcsDepth(end,j))-((ncellsbottom(1,j)-floor(ncellsbottom(1,j)))/2*St.cell(end,j));
                end
               if C_back.zbottomback(uu+1,j)<0.05*V.mcsBed(1,j)
                  C_back.zbottomback(uu+1,j)=nan;
               else
               end
            else
            end
              C_back.zbottomback(uu+2,j)=0.05*V.mcsBed(1,j);
            end
        end
    end
end

 C_back.zbottomback(C_back.zbottomback==0)=nan;
 
%Prediction vector bottom
for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
        C_back.zbottom(1,j)=(V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut2(1,j),j)))/2;%from bottom
        C_back.zpredbottomback(:,j)=linspace(C_back.zbottom(1,j)-1,C_back.zbottom(1,j)+1,20);
    end
end

%Determinate the cellsize of bottom unmeasured zone
for j=1:St.mback
    for uu=1:size(C_back.zbottomback,1)  
        if uu==1
            C_back.cellbottomBack(uu,j)=abs(C_back.zbottomback(uu+1,j)-C_back.zbottomback(uu,j));
        else
            C_back.cellbottomBack(uu,j)=abs(C_back.zbottomback(uu,j)-C_back.zbottomback(uu-1,j));
        end
    end
end

%%
%Position surface points
for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
        ncellssurf(1,j)=(V.mcsDepth(Cut_back.Backcut1(1,j),j)-(St.cell(1,j)/2))/St.cell(1,j);
    end
end

C_back.zsurfback=nan(max(round(ncellssurf)),St.mback);


for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
        
        if V.mcsDepth(Cut_back.Backcut1(1,j),j)<St.cell(1,j)
            C_back.zsurfback(1,j)=V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut1(1,j),j)/2);
        else
            ncellssurf(1,j)=(V.mcsDepth(Cut_back.Backcut1(1,j),j)-(St.cell(1,j)/2))/St.cell(1,j);
            for uu=1:fix(ncellssurf(1,j))            
                C_back.zsurfback(uu,j)=uu*St.cell(1,j)+(V.mcsBed(1,j)-V.mcsDepth(Cut_back.Backcut1(1,j),j));
            end

            if ncellssurf(1,j)-floor(ncellssurf(1,j))~=0
                C_back.zsurfback(uu+1,j)=(V.mcsBed(1,j)-C_back.zsurfback(uu,j))/2+C_back.zsurfback(uu,j);
            else
            end
        end
    else
    end
    C_back.zsurfback(:,j)=flipud(C_back.zsurfback(:,j));
end
    
C_back.zsurfback(C_back.zsurfback==0)=nan;

%prediction vector bottom points
for j=1:St.mback
    if Cut_back.Backcut(1,j)>0
        C_back.zsurf(1,j)=V.mcsBed(1,j)-(V.mcsDepth(Cut_back.Backcut1(1,j),j)/2);
        C_back.zpredsurfback(:,j)=linspace(C_back.zsurf(1,j)-1,C_back.zsurf(1,j)+1,20);
    end
end

%Determinate the cellsize of bottom unmeasured zone

C_back.cellsurfBack=abs(diff(C_back.zsurfback,1));
t=size(C_back.cellsurfBack,1);
if t<2
else
    C_back.cellsurfBack(t,:)=nanmax(C_back.cellsurfBack(1:t-1,:));%
    C_back.cellsurfBack(t+1,:)=nanmax(C_back.cellsurfBack);%add one for the last value
end
end   
