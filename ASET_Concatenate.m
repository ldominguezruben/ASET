function Array=ASET_Concatenate(V,Css,Cut_vel,Cut_back,VelExtra,CssExtra,St,Meas,NoMeas,C_vel,C_back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate all arrays of Velocity, Concentration, Gss and Gw.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Concentration arrays
if isempty(CssExtra)%without extrapolation
    
    Array.Css=Css;
    Array.Depth_Css=V.mcsDepth;

else %with concentration extrapolate
    
    %Determinate Index in different Zones
    for j=1:St.mback
        Array.indexCssbottom(1,j)=nansum(isfinite(CssExtra.surf(:,j)))+Cut_back.Backcut2(1,j)+1;%start bottom zone
        Array.indexCssbottom(2,j)=nansum(isfinite(CssExtra.surf(:,j)))+Cut_back.Backcut2(1,j)+nansum(isfinite(CssExtra.bottom(:,j)));%end bottom zone
        Array.indexCssmeas(1,j)=nansum(isfinite(CssExtra.surf(:,j)))+1;%start measured zone
        Array.indexCssmeas(2,j)=nansum(isfinite(CssExtra.surf(:,j)))+Cut_back.Backcut2(1,j);%end measured zone
        Array.indexCsssurf(1,j)=nansum(isfinite(CssExtra.surf(:,j)));%end measured zone
    end
    
    %Concatenate the Css
     for j=1:St.mback
        t=1;
        l=0;
        for i=1:St.nback+(size(CssExtra.surf,1)-nansum(isfinite(CssExtra.surf(:,j))))
            if i<size(CssExtra.surf,1)+1%Surface zone
                if isnan(CssExtra.surf(i,j))
                    l=l+1;
                else
                    Array.Css(t,j) = CssExtra.surf(i,j);
                    Array.Depth_Css(t,j) = V.mcsBed(1,j) - C_back.zsurfback(i,j);
                    t=t+1;
                end
            elseif i>Array.indexCssmeas(1,j)-1+l & i<Array.indexCssmeas(2,j)+1+l%measured zone
                if isnan(Css(i-(Array.indexCsssurf(1,j)+l),j))
                    %empty
                else
                    Array.Css(t,j) = Css(i-(Array.indexCsssurf(1,j)+l),j);
                    Array.Depth_Css(t,j) = V.mcsBed(1,j) - C_back.zfrombottom(i-(Array.indexCsssurf(1,j)+l),j);
                    t=t+1;
                end
            elseif i> Array.indexCssbottom(1,j)-1+l & i<Array.indexCssbottom(2,j)+1+l%bottom zone
                    if isfinite(CssExtra.bottom(i-(Array.indexCssbottom(1,j)-1+l),j))
                        Array.Css(t,j) = CssExtra.bottom(i-(Array.indexCssbottom(1,j)-1+l),j);
                        Array.Depth_Css(t,j) = V.mcsBed(1,j) - C_back.zbottomback(i-(Array.indexCssbottom(1,j)-1+l),j);
                        t=t+1;
                    else
                        Array.Css(t,j) = nan;
                        Array.Depth_Css(t,j) =  V.mcsBed(1,j) - C_back.zfrombottom(i-(Array.indexCssbottom(1,j)-1+l),j);
                        t=t+1;
                    end
            elseif i> Array.indexCssbottom(2,j)-1+l & i<St.nback+1%rest of zones
                Array.Css(t,j) = nan;
                Array.Depth_Css(t,j) =  V.mcsBed(1,j) - C_back.zfrombottom(i,j);
                t=t+1;
           elseif i> St.nback%rest of zones
                Array.Css(t,j) = nan;
                Array.Depth_Css(t,j) =  nan;
                t=t+1;
            end
        end
        clear t l
     end
    

    %Delete bad data
    Array.Css(Array.Css==Inf)=nan;%Final Array
    
end
    
%%
%Velocity arrays
 if isempty(VelExtra)
     
    Array.uMag=V.mcsMag;%Magnitude Vel
    Array.Depth_uMag=V.mcsDepth;%Depth Vel
    
 else
     
    %Determinate Index differente zones
    for j=1:St.mvel
        Array.indexUbottom(1,j)=nansum(isfinite(VelExtra.surf(:,j)))+Cut_vel.Velcut2(1,j)+1;%start bottom zone
        Array.indexUbottom(2,j)=nansum(isfinite(VelExtra.surf(:,j)))+Cut_vel.Velcut2(1,j)+nansum(isfinite(VelExtra.bottom(:,j)));%end bottom zone
        Array.indexUmeas(1,j)=nansum(isfinite(VelExtra.surf(:,j)))+1;%start measured
        Array.indexUmeas(2,j)=nansum(isfinite(VelExtra.surf(:,j)))+Cut_vel.Velcut2(1,j);%end measured
        Array.indexUsurf(1,j)=nansum(isfinite(VelExtra.surf(:,j)));%end surface zone
    end
    
    
   %Concatenate the Vel
    for j=1:St.mvel
        t=1;
        l=0;
        for i=1:St.nvel+(size(VelExtra.surf,1)-nansum(isfinite(VelExtra.surf(:,j))))
            if i<size(VelExtra.surf,1)+1%Surface zone
                if isnan(VelExtra.surf(i,j))
                    l=l+1;
                else
                    Array.uMag(t,j) = VelExtra.surf(i,j)*100;
                    Array.Depth_uMag(t,j) = V.mcsBed(1,j) - C_vel.zsurfvel(i,j);
                    t=t+1;
                end
            elseif i>Array.indexUmeas(1,j)-1+l & i<Array.indexUmeas(2,j)+1+l%measured zone
                if isnan(V.mcsMag(i-(Array.indexUsurf(1,j)+l),j))
                else
                    Array.uMag(t,j) = V.mcsMag(i-(Array.indexUsurf(1,j)+l),j);
                    Array.Depth_uMag(t,j) = V.mcsBed(1,j) - C_vel.zfrombottom(i-(Array.indexUsurf(1,j)+l),j);
                    t=t+1;
                end
            elseif i> Array.indexUbottom(1,j)-1+l & i<Array.indexUbottom(2,j)+1+l%bottom zone
                    if isfinite(VelExtra.bottom(i-(Array.indexUbottom(1,j)-1+l),j))
                        Array.uMag(t,j) = VelExtra.bottom(i-(Array.indexUbottom(1,j)-1+l),j)*100;
                        Array.Depth_uMag(t,j) = V.mcsBed(1,j) - C_vel.zbottomvel(i-(Array.indexUbottom(1,j)-1+l),j);
                        t=t+1;
                    else
                        Array.uMag(t,j) = nan;
                        Array.Depth_uMag(t,j) =  V.mcsBed(1,j) - C_vel.zfrombottom(i-(Array.indexUbottom(1,j)-1+l),j);
                        t=t+1;
                    end
            elseif i> Array.indexUbottom(2,j)-1+l & i<St.nvel+1%rest of zones
                Array.uMag(t,j) = nan;
                Array.Depth_uMag(t,j) =  V.mcsBed(1,j) - C_vel.zfrombottom(i,j);
                t=t+1;
           elseif i> St.nvel%rest of zones
                Array.uMag(t,j) = nan;
                Array.Depth_uMag(t,j) =  nan;
                t=t+1;
            end
        end
        clear t l
    end
    
    %delete bad data
    Array.uMag(Array.uMag==Inf)=nan;%Final Array Vel Magnitude

 end
 
 %Bed
Array.Bed=V.mcsBed;
%Distance
Array.Dist=V.mcsDist;

%Define size of row
if St.nback<St.nvel
    St.n=St.nback;
else
    St.n=St.nvel;
end

%Define size of column
if St.mback<St.mvel
    St.m=St.mback;
else
    St.m=St.mvel;
end

    
%Concatenate the Qss and Qw array 
 if isempty(VelExtra) | isempty(CssExtra)
     
    Array.Gss = Meas.coarse;
    Array.Gw = Meas.fine;
    
 else
     
    for j=1:St.m
        t=1;
        l=0;
        for i=1:St.n
            if i<size(NoMeas.Surf.coarse,1)+1%Surface zone
                if isnan(NoMeas.Surf.coarse(i,j))
                    l=l+1;
                else
                    Array.Gss(t,j) = NoMeas.Surf.coarse(i,j);
                    Array.Gw(t,j) = NoMeas.Surf.fine(i,j);
                    t=t+1;
                end
            elseif i>Array.indexCssmeas(1,j)-1+l & i<Array.indexCssmeas(2,j)+1+l%measured zone
                if isnan(Meas.coarse(i-(Array.indexCsssurf(1,j)+l),j))
                else
                    Array.Gss(t,j) = Meas.coarse(i-(Array.indexCsssurf(1,j)+l),j);
                    Array.Gw(t,j) = Meas.fine(i-(Array.indexCsssurf(1,j)+l),j);
                    t=t+1;
                end
            elseif i> Array.indexCssbottom(1,j)-1+l & i<Array.indexCssbottom(2,j)+1+l%bottom zone
                    if isfinite(NoMeas.Bottom.coarse(i-(Array.indexCssbottom(1,j)-1+l),j))
                        Array.Gss(t,j) = NoMeas.Bottom.coarse(i-(Array.indexCssbottom(1,j)-1+l),j);
                        Array.Gw(t,j) = NoMeas.Bottom.fine(i-(Array.indexCssbottom(1,j)-1+l),j);
                        t=t+1;
                    else
                        Array.Gss(t,j) = nan;
                        Array.Gw(t,j) = nan;
                        t=t+1;
                    end
            elseif i> Array.indexUbottom(2,j)-1+l%rest of zones
                Array.Gss(t,j) = nan;
                Array.Gw(t,j) = nan;
                t=t+1;
            end
        end
        clear t l
    end
 end

Array.Depth_Gss = Array.Depth_Css;
Array.Depth_Gw = Array.Depth_uMag;

%Graphics results
for j=1:St.mback
    Array.Css_ave(:,j)=Array.Css(:,j)./nanmax(Array.Css(:,j));
    Array.uMag_ave(:,j)=Array.uMag(:,j)./nanmax(Array.uMag(:,j));
    Array.Depth_uMag_ave(:,j)= Array.Depth_uMag(:,j)./Array.Bed(1,j);
    Array.Depth_Css_ave(:,j)= Array.Depth_Css(:,j)./Array.Bed(1,j);
end

%%
% This section find the 20%, 50%, 60% and 80% of the depth and calculate
% the min and max in this range for the Concetration and Velocity

% Velocity
for j=1:St.mvel
    por20=0.2; % 
    por40=0.4; % 
    por60=0.6; % 
    por80=0.8; % 
    
    [~,pos20(j)]=nanmin(abs(Array.Depth_uMag_ave(:,j)-por20));
    [~,pos40(j)]=nanmin(abs(Array.Depth_uMag_ave(:,j)-por40));
    [~,pos60(j)]=nanmin(abs(Array.Depth_uMag_ave(:,j)-por60));
    [~,pos80(j)]=nanmin(abs(Array.Depth_uMag_ave(:,j)-por80));
    
    Array.Vel20(j)=Array.uMag_ave(pos20(j),j);
    Array.Vel40(j)=Array.uMag_ave(pos40(j),j);
    Array.Vel60(j)=Array.uMag_ave(pos60(j),j);
    Array.Vel80(j)=Array.uMag_ave(pos80(j),j);
    
end

%minimum
Array.Velmin=[nanmin(Array.Vel20); nanmin(Array.Vel40); nanmin(Array.Vel60); nanmin(Array.Vel80)];

%maximum
Array.Velmax=[nanmax(Array.Vel20); nanmax(Array.Vel40); nanmax(Array.Vel60);nanmax(Array.Vel80)];

% Concentration
for j=1:St.mback
    por20=0.2; % 
    por40=0.4; % 
    por60=0.6; % 
    por80=0.8; % 
    
    [~,pos20(j)]=min(abs(Array.Depth_Css_ave(:,j)-por20));
    [~,pos40(j)]=min(abs(Array.Depth_Css_ave(:,j)-por40));
    [~,pos60(j)]=min(abs(Array.Depth_Css_ave(:,j)-por60));
    [~,pos80(j)]=min(abs(Array.Depth_Css_ave(:,j)-por80));
    
    Array.Css20(j)=Array.Css_ave(pos20(j),j);
    Array.Css40(j)=Array.Css_ave(pos40(j),j);
    Array.Css60(j)=Array.Css_ave(pos60(j),j);
    Array.Css80(j)=Array.Css_ave(pos80(j),j);
    
end

%minimum
Array.Cssmin=[nanmin(Array.Css20); nanmin(Array.Css40); nanmin(Array.Css60);nanmin(Array.Css80)];

%maximum
Array.Cssmax=[nanmax(Array.Css20); nanmax(Array.Css40); nanmax(Array.Css60);nanmax(Array.Css80)];


