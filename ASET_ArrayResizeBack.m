function [C]=ASET_ArrayResizeBack(V,St)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines the array size of bottom and surface backscatter 
% data. Also determinate the cell on the surface and bottom extrapolation.

% Dominguez Ruben, L. UNL-FICH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:St.mback
    for i=1:St.nback
        C.zfrombottom(i,j)=V.mcsBed(1,j)-V.mcsDepth(i,j);
    end
end

%Last cell with data
%For ascii
V.mcsBack(V.mcsBack==255)=nan;
C.Backcut=sum(isfinite(V.mcsBack));


%Find first cell with data
for j=1:St.mback
    if C.Backcut(1,j)==0
        C.Backcut1(1,j)=0;
    elseif C.Backcut(1,j)==1
        C.Backcut1(1,j)=find(isfinite(V.mcsBack(:,j)));
    elseif C.Backcut(1,j)==2
        findpos=find(isfinite(V.mcsBack(:,j)));
        C.Backcut1(1,j)=min(findpos);
    elseif C.Backcut(1,j)>2
    for i=1:St.nback-2
        if isfinite(V.mcsBack(i,j)) & isfinite(V.mcsBack(i+1,j)) & isfinite(V.mcsBack(i+2,j))
            C.Backcut1(1,j)=i;
            break
        elseif isfinite(V.mcsBack(i,j)) & isfinite(V.mcsBack(i+1,j))
            C.Backcut1(1,j)=i;
            break
        elseif isfinite(V.mcsBack(i,j))
            C.Backcut1(1,j)=i;
            break
        end
    end
    end
end


%Find last cell with data
for j=1:St.mback
    if C.Backcut(1,j)==0%with all data
        C.Backcut2(1,j)=0;
    elseif C.Backcut(1,j)==1%no data on the first cell
        C.Backcut2(1,j)=find(isfinite(V.mcsBack(:,j)));
    elseif C.Backcut(1,j)==2%no data on the second cell
        findpos=find(isfinite(V.mcsBack(:,j)));
        C.Backcut2(1,j)=max(findpos);
    elseif C.Backcut(1,j)>2% the rest of cell
        for i=3:St.nback
         findpos=find(isfinite(V.mcsBack(:,j)));
         last(1,j)=findpos(end);
        if i==last(1,j) & isfinite(V.mcsBack(i,j)) & isfinite(V.mcsBack(i-1,j)) & isfinite(V.mcsBack(i-2,j))
            C.Backcut2(1,j)=i;
            break
        elseif i==last(1,j) & isfinite(V.mcsBack(i,j)) & isfinite(V.mcsBack(i-1,j))
            C.Backcut2(1,j)=i;
            break
        elseif i==last(1,j) & isfinite(V.mcsBack(i,j))
            C.Backcut2(1,j)=i;
            break
        end
        end
    end
end
