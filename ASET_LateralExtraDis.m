function [QEdge]=ASET_LateralExtraDis(Method,Array,V,S,Edge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lateral extrapolation methods. Calculate differente methods of
%extrapolations
%Dominguez Ruben L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes
%Method.Edge 0 triangular
%Method.Edge 1 rectangular


Trifactor=0.3535;%Triangular bank shape see Teledyne (2015)
Recfactor=0.91;%Rectangular bank shape see teledyne (2015)

%find the first and last valid ensemble
%Velocity
for j = 1:size(Array.uMag,2)
    if nansum(isfinite(Array.uMag(:,j)))>0
        Velleft=nanmean(Array.uMag(:,j))/100;
        break   
    else
    end
 
end

for j = size(Array.uMag,2):-1:1
    if nansum(isfinite(Array.uMag(:,j)))>0
        Velright=nanmean(Array.uMag(:,j))/100;
        break 
    end   
end

%Concentration
for j = 1:size(Array.Css,2)
    if nansum(isfinite(Array.Css(:,j)))>0
        Cssleft=nanmean(Array.Css(:,j));
        break 
    end   
end

for j = size(Array.Css,2):-1:1
    if nansum(isfinite(Array.Css(:,j)))>0
        Cssright=nanmean(Array.Css(:,j));
        break  
    end 
end
    
if Method.Edge==0%Triangular
        %Solid discharge
        QEdge.Gss.left=Trifactor*Velleft*Cssleft*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Gss.right=Trifactor*Velright*Cssright*Edge.rightdistance*V.mcsBed(1,end);
        QEdge.Gw.left=Trifactor*Velleft*S.Csf(:,1)/1000*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Gw.right=Trifactor*Velright*S.Csf(:,end)/1000*Edge.rightdistance*V.mcsBed(1,end);
        
        %Liquid Discharge
        QEdge.Q.left=Trifactor*Velleft*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Q.right=Trifactor*Velright*Edge.rightdistance*V.mcsBed(1,end);
        
elseif Method.Edge==1%Rectangular
        %Solid discharge
        QEdge.Gss.left=Recfactor*Velleft*Cssleft*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Gss.right=Recfactor*Velright*Cssright*Edge.rightdistance*V.mcsBed(1,end);
        QEdge.Gw.left=Recfactor*Velleft*S.Csf(:,1)/1000*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Gw.right=Recfactor*Velright/100*S.Csf(:,end)/1000*Edge.rightdistance*V.mcsBed(1,end);
            
        %Liquid Discharge
        QEdge.Q.left=Recfactor*Velleft*Edge.leftdistance*V.mcsBed(1,1);
        QEdge.Q.right=Recfactor*Velright*Edge.rightdistance*V.mcsBed(1,end);
end


%si esta vacio de datos

if isempty(Edge.leftdistance)
    warning('Doesnt exist the left distance to the margin')
    QEdge.Gss.left=0;
    QEdge.Gw.left=0;
    QEdge.Q.left=0;
end
if isempty(Edge.rightdistance)
    warning('Doesnt exist the right distance to the margin')
    QEdge.Gss.right=0;
    QEdge.Gw.right=0;
    QEdge.Q.right=0;
end
