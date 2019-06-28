function Dis=ASET_DischargeLiquid(V,VelExtra,C_vel,St,A,readformatfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Liquid Discharge in the different zones:
% measured, bottom and surface.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(VelExtra)
    if readformatfile==1%VMT
        %Measured zone
        for j=2:size(V.u,2)
            for i=1:size(V.u,1)
                    Dis.Meas(i,j)=(V.mcsMag(i,j)/100)*(V.mcsDist(1,j)-V.mcsDist(1,j-1))*St.cell(i,j);
            end
        end
        
    elseif readformatfile==2%ASCII
        
        Dis.Meas=A.Q.unit;
        
    end
    
    Dis.Total=nansum(nansum(Dis.Meas));
    Dis.bottomT=0;
    Dis.surfT=0;
    
%%
else%with extrapolation methods
    if readformatfile==1%VMT
        %Measured zone
        for j=2:size(V.u,2)
            for i=1:size(V.u,1)
                     Dis.Meas(i,j)=(V.mcsMag(i,j)/100)*(V.mcsDist(1,j)-V.mcsDist(1,j-1))*St.cell(i,j);
            end
        end
    elseif readformatfile==2%ASCII
        Dis.Meas=A.Q.unit;
    end
        
    Dis.MeasT=nansum(nansum(Dis.Meas));
    
    %Unmeasured zone bottom
    for j=2:size(VelExtra.bottom,2)
        for i=1:size(VelExtra.bottom,1)
                Dis.bottom(i,j)=VelExtra.bottom(i,j)*(V.mcsDist(1,j)-V.mcsDist(1,j-1))*C_vel.cellbottomVel(i,j);
        end
    end

    Dis.bottomT=nansum(nansum(Dis.bottom));

    %Unmeasured zone surface
    for j=2:size(VelExtra.surf,2)
        for i=1:size(VelExtra.surf,1)
                Dis.surf(i,j)=VelExtra.surf(i,j)*(V.mcsDist(1,j)-V.mcsDist(1,j-1))*C_vel.cellsurfVel(i,j);
        end
    end

    Dis.surfT=nansum(nansum(Dis.surf));

    Dis.Total=Dis.MeasT+Dis.bottomT+Dis.surfT;

end

     

  


  