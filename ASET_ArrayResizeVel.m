function [C]=ASET_ArrayResizeVel(V,St)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function defines the array size of bottom and surface velocity data.

%Dominguez Ruben,L. UNL-FICH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:St.mvel
    for i=1:St.nvel
        C.zfrombottom(i,j)=V.mcsBed(1,j)-V.mcsDepth(i,j);
    end
end

%Cut variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last cell with data
C.Velcut=nansum(isfinite(V.mcsMag));

%First cell with data
for j=1:St.mvel   
    if C.Velcut(1,j)==0
        C.Velcut1(1,j)=0;
    elseif C.Velcut(1,j)==1
        C.Velcut1(1,j)=find(isfinite(V.mcsMag(:,j)));
    elseif C.Velcut(1,j)==2
        findpos=find(isfinite(V.mcsMag(:,j)));
        C.Velcut1(1,j)=min(findpos);
    elseif C.Velcut(1,j)>2
        for i=1:St.nvel-2
            if isfinite(V.mcsMag(i,j)) & isfinite(V.mcsMag(i+1,j)) & isfinite(V.mcsMag(i+2,j))
                C.Velcut1(1,j)=i;
            break
            elseif isfinite(V.mcsMag(i,j)) & isfinite(V.mcsMag(i+1,j))
                C.Velcut1(1,j)=i;
            break
            elseif isfinite(V.mcsMag(i,j))
                C.Velcut1(1,j)=i;
            break
            end
        end
    end
end


%last cell with data
for j=1:St.mvel
    if C.Velcut(1,j)==0
        C.Velcut2(1,j)=0;
    elseif C.Velcut(1,j)==1
        C.Velcut2(1,j)=find(isfinite(V.mcsMag(:,j)));
    elseif C.Velcut(1,j)==2
        findpos=find(isfinite(V.mcsMag(:,j)));
        C.Velcut2(1,j)=max(findpos);
    elseif C.Velcut(1,j)>2
        for i=3:St.nvel
         findpos=find(isfinite(V.mcsMag(:,j)));
         last(1,j)=findpos(end);
        if i==last(1,j) & isfinite(V.mcsMag(i,j)) & isfinite(V.mcsMag(i-1,j)) & isfinite(V.mcsMag(i-2,j))
            C.Velcut2(1,j)=i;
            break
        elseif i==last(1,j) & isfinite(V.mcsMag(i,j)) & isfinite(V.mcsMag(i-1,j))
            C.Velcut2(1,j)=i;
            break
        elseif i==last(1,j) & isfinite(V.mcsMag(i,j))
            C.Velcut2(1,j)=i;
            break
        end
        end
    end
end