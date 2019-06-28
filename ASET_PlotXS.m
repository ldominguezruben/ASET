function [fig]=ASET_PlotXS(cof,Method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plot the Cross section contour map. The variables possible
% are Velocity, Ms2, Gss and Gw.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot Cross-Section

%Load data
guiparams = getappdata(0,'guiparams');

VelExtraM= guiparams.VelExtraM;
CssExtraM=guiparams.CssExtraM;

 V   = guiparams.V;
Meas   = guiparams.Meas;
Css   = guiparams.Css;
Array   = guiparams.Array;


USys = getappdata(0,'units');

switch cof.crosssel
    case 1 %Css
        
        Array.Css=smooth2a(Array.Css,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_Css=smooth2a(Array.Depth_Css,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(Array.Dist,cof.horizontalsmooth,cof.verticalsmooth);

        f=figure(1); 
        fig.position=set(f, 'Position', [10 10 1000 500]);
        pcolor(Array.Dist,Array.Depth_Css,Array.Css*1000);
        shading interp
        if USys==0
        title('ASET: Mean Cross Section Contour [mg/l]','Fontsize',18)
        xlabel('Distance across section [m]','Fontsize',16);ylabel('Depth [m]','Fontsize',16);
        fig.Label.String = 'Ms2 [mg/l]';
        elseif USys==1
        title('ASET: Mean Cross Section Contour [t/ft3]','Fontsize',18)
        xlabel('Distance across section [ft]','Fontsize',16);ylabel('Depth [ft]','Fontsize',16);
        fig.Label.String = 'Ms2 [t/ft3]';
        end
        
        grid off
        hold on
        plot(Array.Dist(1,1:end),Array.Bed,'-k','Linewidth',2)
        hold off
        set(gca,'ydir','reverse');
        if Method==1;
            ylim([0 max(max(Array.Bed))]);
            fig.xlim=get(gca,'xlim');
            fig.ylim=get(gca,'ylim');
            fig.zmin=1.35*floor(nanmin(nanmin(Array.Css*1000)));
            fig.zmax=0.25*ceil(nanmax(nanmax(Array.Css*1000)));
            caxis([fig.zmin fig.zmax]);
            fig.colormap=colormap('jet');
            fig.colorbar=colorbar;
        else
            fig.xlim=set(gca,'xlim',[cof.minxaxis cof.maxxaxis]);
            fig.ylim=set(gca,'ylim',[cof.minyaxis cof.maxyaxis]);
            caxis([cof.minvel cof.maxvel]);
            fig.colormap=colormap(cof.colormaps_mcs);
            fig.colorbar=colorbar;
        end
        fig.FontSize=16;


    case 2 %GSS

      %  Array.Gss(Array.Gss==0)=nan;
        Array.Gss=smooth2a(Array.Gss,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_Gss=smooth2a(Array.Depth_Gss,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(Array.Dist,cof.horizontalsmooth,cof.verticalsmooth);

        f=figure(1); 
        fig.position=set(f, 'Position', [10 10 1000 500]);
        pcolor(Array.Dist,Array.Depth_Gss,Array.Gss);shading interp
        
        if USys==0
            title('ASET: Mean Cross Section Contour [kg/s]','Fontsize',18)
            xlabel('Cross Section [m]','Fontsize',16);ylabel('Depth [m]','Fontsize',16)
            fig.Label.String = 'Suspended Bed Sediment Transport Qss[kg/s]';
        elseif USys==1
            title('ASET: Mean Cross Section Contour [t/s]','Fontsize',18)
            xlabel('Cross Section [ft]','Fontsize',16);ylabel('Depth [ft]','Fontsize',16)
            fig.Label.String = 'Suspended Bed Sediment Transport Qss[t/s]';
        end
        
        grid off
        hold on
        plot(Array.Dist,Array.Bed,'-k','Linewidth',2)
        fig.ydir=set(gca,'ydir','reverse');
        hold off
       if Method==1;
            ylim([0 max(max(Array.Bed))]);
            fig.xlim=get(gca,'xlim');
            fig.ylim=get(gca,'ylim');
            fig.zmin=0;%
            fig.zmax=nanmax(nanmax(Array.Gss));
            caxis([fig.zmin fig.zmax])
            fig.colormap=colormap('jet');
            fig.colorbar=colorbar;
        else
            fig.xlim=set(gca,'xlim',[cof.minxaxis cof.maxxaxis]);
            fig.ylim=set(gca,'ylim',[cof.minyaxis cof.maxyaxis]);
            caxis([cof.minvel cof.maxvel]);
            fig.colormap=colormap(cof.colormaps_mcs);
            fig.colorbar=colorbar;
       end
        
        fig.FontSize=12;
    case 3 %Gw

        Meas.fine(Meas.fine==0)=nan;
        Array.Gw=smooth2a(Array.Gw,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_Gw=smooth2a(Array.Depth_Gw,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(Array.Dist,cof.horizontalsmooth,cof.verticalsmooth);


        f=figure(1); 
        fig.position=set(f, 'Position', [10 10 1000 500]);
        pcolor(Array.Dist,Array.Depth_Gw,Array.Gw);shading interp
        
        if USys==0;
            title('ASET: Mean Cross Section Contour [kg/s]','Fontsize',18)
            xlabel('Cross Section [m]','Fontsize',16);ylabel('Depth [m]','Fontsize',16)
            fig.Label.String = 'Washload Qw [kg/s]';
        elseif USys==1
            title('ASET: Mean Cross Section Contour [t/s]','Fontsize',18)
            xlabel('Cross Section [ft]','Fontsize',16);ylabel('Depth [ft]','Fontsize',16)
            fig.Label.String = 'Washload Qw [t/s]';
        end
        
        grid off
        hold on
        plot(Array.Dist,Array.Bed,'-k','Linewidth',2)
        hold off
        fig.ydir=set(gca,'ydir','reverse');
        
       if Method==1;
            ylim([0 max(max(Array.Bed))]);
            fig.xlim=get(gca,'xlim');
            fig.ylim=get(gca,'ylim');
            fig.zmin=nanmin(nanmin(Array.Gw));
            fig.zmax=0.952*ceil(nanmax(nanmax(Array.Gw)));
            caxis([fig.zmin fig.zmax])
            fig.colormap=colormap('jet');
            fig.colorbar=colorbar;
        else
            fig.xlim=set(gca,'xlim',[cof.minxaxis cof.maxxaxis]);
            fig.ylim=set(gca,'ylim',[cof.minyaxis cof.maxyaxis]);
            caxis([cof.minvel cof.maxvel]);
            fig.colormap=colormap(cof.colormaps_mcs);
            fig.colorbar=colorbar;
       end
        
        fig.FontSize=12;  
    case 4 %Velocity

        Array.uMag=smooth2a(Array.uMag,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_uMag=smooth2a(Array.Depth_uMag,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(Array.Dist,cof.horizontalsmooth,cof.verticalsmooth);

        f=figure(1); 
        fig.position=set(f, 'Position', [10 10 1000 500]);
        pcolor(Array.Dist,Array.Depth_uMag,Array.uMag);shading interp
        
        if USys==0
            title('ASET: Mean Cross Section Contour [cm/s]','Fontsize',18)
            xlabel('Cross Section [m]','Fontsize',16);ylabel('Depth [m]','Fontsize',16)
            fig.Label.String = 'Streamwise Velocity [cm/s]';
        elseif USys==1
            title('ASET: Mean Cross Section Contour [ft/s]','Fontsize',18)
            xlabel('Cross Section [ft]','Fontsize',16);ylabel('Depth [ft]','Fontsize',16)
            fig.Label.String = 'Streamwise Velocity [ft/s]';
        end
        
        grid off
        hold on
        plot(Array.Dist,Array.Bed,'-k','Linewidth',2)
        fig.ydir=set(gca,'ydir','reverse');
        if Method==1;%first time
            ylim([0 max(max(Array.Bed))]);
            fig.xlim=get(gca,'xlim');
            fig.ylim=get(gca,'ylim');
            fig.zmin=0;%
            fig.zmax=ceil(nanmax(nanmax(Array.uMag)));
            caxis([fig.zmin fig.zmax])
            fig.colormap=colormap('jet');
            fig.colorbar=colorbar;
        else
            fig.xlim=set(gca,'xlim',[cof.minxaxis cof.maxxaxis]);
            fig.ylim=set(gca,'ylim',[cof.minyaxis cof.maxyaxis]);
            caxis([cof.minvel cof.maxvel]);
            fig.colormap=colormap(cof.colormaps_mcs);
            fig.colorbar=colorbar;
        end
        
        fig.FontSize=12;
end

if Method==1
    fig.exaggeration=floor(max(fig.xlim)/max(fig.ylim));
    if fig.exaggeration==0
        fig.exaggeration=1;
    else
    end
    
    fig.smoothvertical=0;
    fig.smoothhorizontal=0;
else
    fig.exaggeration=cof.verticalexaggeration;
    fig.smoothvertical=cof.verticalsmooth;
    fig.smoothhorizontal=cof.horizontalsmooth;
end

% Adjust the plot
set(gca,...
    'DataAspectRatio',   [fig.exaggeration 1 1],...
    'PlotBoxAspectRatio',[fig.exaggeration 1 1]...
    ...'FontSize',          14)
    )

