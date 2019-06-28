function [VelExtra]=ASET_ExtrapolVel(V,Cut_vel,St,INT,C_vel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the extrapolation value for different methods in
% velocity data: Constant value, 3pt data and Law of the Wall.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surfsize=size(C_vel.zsurfvel,1);
bottomsize=size(C_vel.zbottomvel,1);

for j=1:St.mvel
    if Cut_vel.Velcut(1,j)<=6%default neccesary data to extrapolated  
                                
        if Cut_vel.Velcut1(1,j)==0 
            VelExtra.surf(1:surfsize,j)=nan;
        else   
            for e=1:length(C_vel.zsurfvel(:,j));
                if isnan(C_vel.zsurfvel(e,j))
                    VelExtra.surf(e,j)=nan;
                else
                    VelExtra.surf(e,j)=V.mcsMag(Cut_vel.Velcut1(1,j),j)/100; %extrapolate constant first valid data    
                end
            end

         end

        if Cut_vel.Velcut2(1,j)==0 
            VelExtra.bottom(1:bottomsize,j)=nan;
        else
            for e=1:length(C_vel.zbottomvel(:,j));
                if isnan(C_vel.zbottomvel(e,j))
                    VelExtra.bottom(e,j)=nan;
                else
                    VelExtra.bottom(e,j)=V.mcsMag(Cut_vel.Velcut2(1,j),j)/100;% extrapolated constat last data
                end
            end       
        end
            
    else %extrapolate methods            
            if (strcmp(INT,'CON'))% 
                %%%%%%%%%%%%%%%%%%%%%%%
                %Constant extrapolation
                %%%%%%%%%%%%%%%%%%%%%%%
                [VelEx]=ASET_ExtrapVConstant(V,C_vel,j,Cut_vel);

                VelExtra.surf(:,j) = VelEx.surf;

                VelExtra.bottom(:,j) =  VelEx.bottom;

            elseif (strcmp(INT,'LVI'))%

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %Linear Velocity extrapolation 3pt
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [VelEx]=ASET_ExtrapV3pt(V,C_vel,j,Cut_vel);

                VelExtra.surf(:,j) = VelEx.surf;
                VelExtra.predsurf(:,j)=VelEx.predsurf;

                VelExtra.bottom(:,j)=VelEx.bottom;
                VelExtra.predbottom(:,j)=VelEx.predbottom;

            elseif (strcmp(INT,'LPVI'))

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %No-Linear Velocity Extrapolation Law of the Wall
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               [VelEx]=ASET_ExtrapVLawOfTheWall(V,C_vel,j,Cut_vel);

               VelExtra.pred(1:length(VelEx.zpred),j) = VelEx.pred;
               VelExtra.zpred(:,j) = VelEx.zpred;
               VelExtra.surf(1:surfsize,j) = VelEx.surf;
               VelExtra.bottom(1:bottomsize,j)= VelEx.bottom;
               VelExtra.ustarwl(j)=VelEx.ustarwl;
               VelExtra.ordenada(j)=VelEx.ordenada;
               VelExtra.rcuad(j)=VelEx.rsqvel;

            end

            %Condition if extrapolation create image numbers
            if isreal(VelExtra.surf(1:surfsize,j)) |  isreal(VelExtra.bottom(1:bottomsize,j))
            else %return the loop to constant regresion
                
                [VelEx]=ASET_ExtrapVConstant(V,C_vel,j,Cut_vel);

                 VelExtra.surf(:,j) = VelEx.surf;
                 VelExtra.bottom(:,j) = VelEx.bottom;
            end
    end
                        
end

%Control the calculate data
VelExtra.surf(VelExtra.surf==Inf)=nan;
VelExtra.surf(VelExtra.surf==-Inf)=nan;
VelExtra.bottom(VelExtra.bottom==Inf)=nan;    
VelExtra.bottom(VelExtra.bottom==-Inf)=nan;    