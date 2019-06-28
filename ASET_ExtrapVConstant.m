function [VelExtra]=ASET_ExtrapVConstant(V,C_vel,j,Cut_vel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extrapolate the velocity variable using the first and last
% value of the ensemble.
%
% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Surface zone velocity
for e=1:length(C_vel.zsurfvel(:,j));
    if isnan(C_vel.zsurfvel(e,j))
        VelExtra.surf(e,1)=nan;
    else
        VelExtra.surf(e,1)=V.mcsMag(Cut_vel.Velcut1(1,j),j)/100;       
    end
end

%Bottom zone velocity
for e=1:length(C_vel.zbottomvel(:,j));
    if isnan(C_vel.zbottomvel(e,j))
        VelExtra.bottom(e,1)=nan;
    else
        VelExtra.bottom(e,1)=V.mcsMag(Cut_vel.Velcut2(1,j),j)/100;
    end
end   
    
    
    