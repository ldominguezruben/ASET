function [Meas]=ASET_TransportMeas(Css,V,S,St,Dis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the value of Qss in the measured zone.

%by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check the min array size

if St.nback<St.nvel
    St.n=St.nback;
else
    St.n=St.nvel;
end

if St.mback<St.mvel
    St.m=St.mback;
else
    St.m=St.mvel;
end

%%
%%%%%%%%%%%%%%%%%%%%%
%Transport calculate
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Qss measured zone Kg/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for j=1:St.m
    for i=1:St.n
     if isnan(V.mcsMag(i,j)) | isnan(V.mcsBack(i,j))
        Meas.coarse(i,j)=nan;%zero is is nan velocity or backsccatter
     else
        Meas.coarse(i,j)=(Css(i,j))*Dis.Meas(i,j);%kg/s
     end  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Qw measured zone kg/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:St.m
    for i=1:St.n
     if isnan(V.mcsMag(i,j)) | isnan(V.mcsBack(i,j))
        Meas.fine(i,j)=0;%zero is is nan velocity or backsccatter
     else
        Meas.fine(i,j)=(S.Csf/1000)*Dis.Meas(i,j);%kg/s
     end
    end
end


