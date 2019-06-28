function varargout = ASET_PlotTest(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Test is a function that verified the extrapolation method aplied. 
% Show the results in the different zone of calculate. Also show it the 
% ustar, Rouse Number and rcuad.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_PlotTest_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_PlotTest_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ASET_PlotTest is made visible.
function ASET_PlotTest_OpeningFcn(hObject, eventdata, handles, varargin)
%Initial Panel
set(handles.figure1,'Name',['ASET: Check Extrapolation '], ...
    'DockControls','off')

%Load data
guiparams=getappdata(0, 'guiparams');
St   = guiparams.St;
Array   = guiparams.Array;
    
%Units 
USys = getappdata(0,'units');

%Value of slider
set(handles.sedimentslider,'Max',St.mback,'Min',0)
svsed=fix(get(handles.sedimentslider,'Value'));

set(handles.velocityslider,'Max',St.mvel,'Min',0)
svvel=fix(get(handles.velocityslider,'Value'));

set(handles.popupVel,'Enable','off');
set(handles.popupCss,'Enable','off');

%[1] Concentraton plot
axes(handles.axesconc)
pcolor(Array.Dist,Array.Depth_Css,Array.Css*1000);shading interp
grid off
hold on
plot(Array.Dist,Array.Bed,'-k','Linewidth',1)
if svsed==0
    line([Array.Dist(1,1) Array.Dist(1,1)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
else
    line([Array.Dist(1,svsed) Array.Dist(1,svsed)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
end
set(gca,'ydir','reverse')
ylim([0 max(max(Array.Bed))])
zmin=1.35*floor(nanmin(nanmin(Array.Css*1000)));
zmax=0.25*ceil(nanmax(nanmax(Array.Css*1000)));
caxis([zmin zmax])
            
if USys==0;
    title('Ms2 [mg/l]','Fontsize',12)
    ylabel('Depth [m]','Fontsize',10)
elseif USys==1;
    title('Ms2 [t/ft3]','Fontsize',12)
    ylabel('Depth [ft]','Fontsize',10)
end

colormap('jet');
cc=colorbar;
cc.FontSize=10;

%[2] Velocity plot
axes(handles.axesvel)
%figure 
        
pcolor(Array.Dist,Array.Depth_uMag,real(Array.uMag./100));shading interp
grid off
hold on
plot(Array.Dist,Array.Bed,'-k','Linewidth',1)
if svvel==0
    line([Array.Dist(1,1) Array.Dist(1,1)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
else
    line([Array.Dist(1,svvel) Array.Dist(1,svvel)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
end
set(gca,'ydir','reverse')
ylim([0 max(max(Array.Bed))])
zmin=floor(nanmin(nanmin(real(Array.uMag))))./100;
zmax=ceil(nanmax(nanmax(real(Array.uMag))))./100;
caxis([zmin zmax])

if USys==0;
    title('Velocity Magnitude [m/s]','Fontsize',12)
    ylabel('Depth [m]','Fontsize',10)
elseif USys==1;
    title('Velocity Magnitude [ft/s]','Fontsize',12)
    ylabel('Depth [ft]','Fontsize',10)
end
colormap('jet');
cc=colorbar;
cc.FontSize=10;
    
 %Show Initial Distance (0)
if USys==0;%Metric units
    if svvel==0
        distance = uicontrol('style','text',...
      'unit','pix',...
      'position',[540 448 150 20],...
      'fontsize',10);
        set(distance,'string',['Distance [m]= 0']);
        drawnow;
        distance = uicontrol('style','text',...
      'unit','pix',...
      'position',[120 448 150 20],...
      'fontsize',10);
        set(distance,'string',['Distance [m]= 0']);
        drawnow; 
    end
elseif USys==1%English units
    if svvel==0
        distance = uicontrol('style','text',...
      'unit','pix',...
      'position',[540 448 150 20],...
      'fontsize',10);
        set(distance,'string',['Distance [ft]= 0']);
        drawnow;
        distance = uicontrol('style','text',...
      'unit','pix',...
      'position',[120 448 150 20],...
      'fontsize',10);
        set(distance,'string',['Distance [ft]= 0']);
        drawnow;
    end
end

axes(handles.axesconcp)
grid on
axes(handles.axesvelp)
grid on

% Choose default command line output for ASET_PlotTest
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ASET_PlotTest_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%Toolbar Menu
%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function filefuction_Callback(hObject, eventdata, handles)
%Empty


% --------------------------------------------------------------------
function exitfunction_Callback(hObject, eventdata, handles)
close

%%%%%%%%%%%%%%
%Main Panel
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Concetration
 
% --- Executes during object creation, after setting all properties.
function sedimentslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sedimentslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function r2conc_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function r2conc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupCss.
function popupCss_Callback(hObject, eventdata, handles)
sedimentslider_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupCss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupCss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sedimentslider_Callback(hObject, eventdata, handles)

%Load data
guiparams=getappdata(0, 'guiparams');
St   = guiparams.St;
CssExtra   = guiparams.CssExtra;
Array   = guiparams.Array;
CssExtraM   = guiparams.CssExtraM;
C_back   = guiparams.C_back;

%Units 
USys = getappdata(0,'units');

%Slider data
set(handles.sedimentslider,'Max',St.mback,'Min',0)
svsed=round(get(handles.sedimentslider,'Value'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distance slider Sed
%%%%%%%%%%%%%%%%%%%%%%%%%%
if USys==0;
    if svsed==0
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[120 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [m]=',num2str(round(Array.Dist(1,1)*100)/100)]);
        drawnow;
    else
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[120 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [m]=',num2str(round(Array.Dist(1,svsed)*100)/100)]);
        drawnow;
    end
elseif USys==1
    if svsed==0
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[120 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [ft]=',num2str(round(Array.Dist(1,1)*100)/100)]);
        drawnow;
    else
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[120 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [ft]=',num2str(round(Array.Dist(1,svsed)*100)/100)]);
        drawnow;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[1] Concentration contour color plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
axes(handles.axesconc)
pcolor(Array.Dist,Array.Depth_Css,Array.Css*1000);shading interp
grid off
hold on
plot(Array.Dist,Array.Bed,'-k','Linewidth',1)
if svsed==0
    line([Array.Dist(1,1) Array.Dist(1,1)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
else
    line([Array.Dist(1,svsed) Array.Dist(1,svsed)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
end
zmin=1.35*floor(nanmin(nanmin(Array.Css*1000)));
zmax=0.25*ceil(nanmax(nanmax(Array.Css*1000)));
caxis([zmin zmax])
set(gca,'ydir','reverse')
ylim([0 max(max(Array.Bed))])
hold off
    
if USys==0
    title('Ms2 [mg/l]','Fontsize',12)
    ylabel('Depth [m]','Fontsize',10)
elseif USys==1
    title('Ms2 [t/ft3]','Fontsize',12)
    ylabel('Depth [ft]','Fontsize',10)
end
colormap('jet');
cc=colorbar;
cc.FontSize=10;

%the first case
if svsed==0
   svsed=1; 
else 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Concentration Profile
%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(CssExtraM,'NONE'))
        
    set(handles.popupCss,'Value',2);
    set(handles.popupCss,'Enable','off');
    
    axes(handles.axesconcp)
    plot(Array.Css(:,svsed)./nanmax(Array.Css(:,svsed)),...
     Array.Depth_Css(:,svsed)./Array.Bed(1,svsed),'ok')
    hold on     
    set(handles.axesconcp,'ydir','reverse');
    title('Ms2 Profile','Fontsize',14)
    xlabel('Ms2/Ms2max','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2conc,'Enable','off');
    set(handles.r2conc,'String','');

elseif (strcmp(CssExtraM,'CON'))
    % Constant extrapolation method
    axes(handles.axesconcp)
    set(handles.popupCss,'Enable','off');
    plot(Array.Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
     Array.Depth_Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)./Array.Bed(1,svsed),'ok')
    hold on     
    set(handles.axesconcp,'ydir','reverse');

    plot(Array.Css(1:Array.indexCsssurf(1,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
     Array.Depth_Css(1:Array.indexCsssurf(1,svsed),svsed)./Array.Bed(1,svsed),'oy','MarkerSize',9,'MarkerFaceColor','y')

    plot(Array.Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
     Array.Depth_Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)./Array.Bed(1,svsed),'ob','MarkerSize',9,'MarkerFaceColor','b')
    
    title('Ms2 Profile','Fontsize',14)
    xlabel('Ms2/Ms2max','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone','Surface Extrapolation','Bottom Extrapolation');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2conc,'Enable','off');
    set(handles.r2conc,'String','');
    
elseif (strcmp(CssExtraM,'LSI'))
    % 3Pt. Slope extrapolation method
    set(handles.popupCss,'Value',2);
    set(handles.popupCss,'Enable','off');
    
    axes(handles.axesconcp)
    
     plot(Array.Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
         Array.Depth_Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)./Array.Bed(1,svsed),'ok')
     hold on     
     set(handles.axesconcp,'ydir','reverse');
     
     plot(CssExtra.predsurf(:,svsed)./nanmax(Array.Css(:,svsed)),...
         (Array.Bed(1,svsed)-C_back.zpredsurfback(:,svsed))./Array.Bed(1,svsed),'-g')
     
     plot(CssExtra.predbottom(:,svsed)./nanmax(Array.Css(:,svsed)),...
         (Array.Bed(1,svsed)-C_back.zpredbottomback(:,svsed))./Array.Bed(1,svsed),'-g')

     plot(Array.Css(1:Array.indexCsssurf(1,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
         Array.Depth_Css(1:Array.indexCsssurf(1,svsed),svsed)./Array.Bed(1,svsed),'oy','MarkerSize',9,'MarkerFaceColor','y')

     plot(Array.Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)./nanmax(Array.Css(:,svsed)),...
         Array.Depth_Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)./Array.Bed(1,svsed),'ob','MarkerSize',9,'MarkerFaceColor','b')
    
    title('Ms2 Profile','Fontsize',14)
    xlabel('Ms2/Ms2max','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone','Fitting Surf','Fitting Bot','Surface Extrapolation','Bottom Extrapolation');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2conc,'Enable','off');
    set(handles.r2conc,'String','');

elseif (strcmp(CssExtraM,'RSI'))%Rouse distribution
         
    axes(handles.axesconcp)
    set(handles.popupCss,'Enable','on');
    sel=get(handles.popupCss,'Value');

    switch sel
        case 1% log vs log plot
           if svsed==0
                cla(handles.axesconcp)
           else
           end

            loglog((Array.Depth_Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)*CssExtra.a(1,svsed))./((Array.Bed(1,svsed)-CssExtra.a(1,svsed))*(Array.Bed(1,svsed)-Array.Depth_Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed))),...
                Array.Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed).*1000,'ok')
            hold on     

            if USys==0
                loglog(CssExtra.acoezpred(:,svsed)',flipud(CssExtra.pred(:,svsed)').*1000,'g')
            elseif USys==1
                loglog(CssExtra.acoezpred(:,svsed)'/0.3048,flipud(CssExtra.pred(:,svsed)'/(907200/28.32)),'g')
            end

            loglog((Array.Depth_Css(1:Array.indexCsssurf(1,svsed),svsed)*CssExtra.a(1,svsed))./((Array.Bed(1,svsed)-CssExtra.a(1,svsed))*(Array.Bed(1,svsed)-Array.Depth_Css(1:Array.indexCsssurf(1,svsed),svsed))),...
                Array.Css(1:Array.indexCsssurf(1,svsed),svsed).*1000,'oy','MarkerSize',9,'MarkerFaceColor','y')

            loglog((Array.Depth_Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)*CssExtra.a(1,svsed))./((Array.Bed(1,svsed)-CssExtra.a(1,svsed))*(Array.Bed(1,svsed)-Array.Depth_Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed))),...
                Array.Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed).*1000,'ob','MarkerSize',9,'MarkerFaceColor','b')

            if USys==0
                loglog(1,CssExtra.Ca(1,svsed).*1000,'or','MarkerSize',9,'MarkerFaceColor','r')
            elseif USys==1
                loglog(((Array.Bed(1,svsed)-CssExtra.a(1,svsed)/(907200/28.32))./CssExtra.a(1,svsed)/(907200/28.32))./((Array.Bed(1,svsed)-CssExtra.a(1,svsed)/(907200/28.32))/CssExtra.a(1,svsed)/(907200/28.32)),...
                    CssExtra.Ca(1,svsed)/(907200/28.32),'or','MarkerSize',9,'MarkerFaceColor','r')
            end

            hold off          
            title('Ms2 Profile','Fontsize',14)
            xlabel('log((H-z)zref/(H-zref)z)','Fontsize',12);ylabel('log(Ms2)','Fontsize',12)
            leg=legend('Measured Zone','Fitting','Surface Extrapolation','Bottom Extrapolation','Cref');
            set(leg,'FontSize',9,'Location','best');
            set(handles.r2conc,'Enable','on');
                %write rcuad
            if CssExtra.rcuad(1,svsed)<0.7
                set(handles.r2conc,'String',num2str(round(CssExtra.rcuad(1,svsed)*100)/100));
                set(handles.r2conc,'BackgroundColor','red');
            else
                set(handles.r2conc,'String',num2str(round(CssExtra.rcuad(1,svsed)*100)/100));
                set(handles.r2conc,'BackgroundColor','white');
            end

        case 2% adimensional plot
           if svsed==0
                cla(handles.axesconcp)
           else
           end

            plot(Array.Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed)./CssExtra.Ca(1,svsed),...
                (Array.Depth_Css(Array.indexCssmeas(1,svsed):Array.indexCssmeas(2,svsed),svsed))./(Array.Bed(1,svsed)),'ok')
            hold on  
            set(handles.axesconcp,'ydir','reverse');
            plot(fliplr(CssExtra.pred(:,svsed)')./(CssExtra.Ca(1,svsed)),...
                CssExtra.zpred(:,svsed)'./Array.Bed(1,svsed),'g');
            plot(Array.Css(1:Array.indexCsssurf(1,svsed),svsed)./CssExtra.Ca(1,svsed),...
                (Array.Depth_Css(1:Array.indexCsssurf(1,svsed),svsed))./(Array.Bed(1,svsed)),'oy','MarkerSize',9,'MarkerFaceColor','y')
            plot(Array.Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed)./CssExtra.Ca(1,svsed),...
                (Array.Depth_Css(Array.indexCssbottom(1,svsed):Array.indexCssbottom(2,svsed),svsed))./(Array.Bed(1,svsed)),'ob','MarkerSize',9,'MarkerFaceColor','b')
            plot(1,0.95,'or','MarkerSize',9,'MarkerFaceColor','r')
            hold off
            title('Ms2 Profile','Fontsize',14)
            xlabel('Ms2/Cref','Fontsize',12);ylabel('z/H','Fontsize',12)
            leg=legend('Measured Zone','Fitting','Surface Extrapolation','Bottom Extrapolation','Cref');
            set(leg,'FontSize',9,'Location','best');
            set(handles.r2conc,'Enable','off');
            set(handles.r2conc,'String','');

    end
            
    guidata(hObject,handles)

    

end
    


% ---------------------
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Velocity

% --- Executes during object creation, after setting all properties.
function velocityslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velocityslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function r2vel_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function r2vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupVel.
function popupVel_Callback(hObject, eventdata, handles)
velocityslider_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function popupVel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function velocityslider_Callback(hObject, eventdata, handles)

%Load data
guiparams=getappdata(0, 'guiparams');
St   = guiparams.St;
VelExtra   = guiparams.VelExtra;
Array   = guiparams.Array;
VelExtraM   = guiparams.VelExtraM;
C_vel = guiparams.C_vel;

%Units 
USys = getappdata(0,'units');


set(handles.velocityslider,'Max',St.mvel,'Min',0)
svvel=fix(get(handles.velocityslider,'Value'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[1] Velocity contour color plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axesvel)
%figure 
pcolor(Array.Dist,Array.Depth_uMag,real(Array.uMag./100));shading interp
grid off
hold on
plot(Array.Dist,Array.Bed,'-k','Linewidth',1)
if svvel==0
    line([Array.Dist(1,1) Array.Dist(1,1)] ,[0 max(max(Array.Bed))],'color','r','LineWidth',2)
else
    line([Array.Dist(1,svvel) Array.Dist(1,svvel)] ,[0 nanmax(nanmax(Array.Bed))],'color','r','LineWidth',2)
end
set(gca,'ydir','reverse')
ylim([0 nanmax(nanmax(Array.Bed))])
zmin=floor(nanmin(nanmin(real(Array.uMag))))./100;
zmax=ceil(nanmax(nanmax(real(Array.uMag))))./100;
caxis([zmin zmax])
    
if USys==0;
    title('Velocity Magnitude [m/s]','Fontsize',12)
    ylabel('Depth [m]','Fontsize',10)
elseif USys==1;
    title('Velocity Magnitude [ft/s]','Fontsize',12)
    ylabel('Depth [ft]','Fontsize',10)
end

hold off
        
colormap('jet');
cc=colorbar;
cc.FontSize=10;
    
%%%%%%%%%%%%%%%%%%%
%Distance Vel
%%%%%%%%%%%%%%%%%%%
if USys==0;
    if svvel==0
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[540 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [m]=',num2str(round(Array.Dist(1,1)*100)/100)]);
        drawnow;
    else
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[540 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [m]=',num2str(round(Array.Dist(1,svvel)*100)/100)]);
        drawnow;
    end
elseif USys==1
    if svvel==0
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[540 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [ft]=',num2str(round(Array.Dist(1,1)*100)/100)]);
        drawnow;
    else
        distance = uicontrol('style','text',...
          'unit','pix',...
          'position',[540 448 150 20],...
          'fontsize',10);
        set(distance,'string',['Distance [ft]=',num2str(round(Array.Dist(1,svvel)*100)/100)]);
        drawnow;
    end
end


if svvel==0    
    svvel=1;
else
end
    
%%%%%%%%%%%%%%%%%%%%%%%%
%Profile Velocity Graph
%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(VelExtraM,'NONE'))
    %None extrapolation Method
    
    set(handles.popupVel,'Value',2);
    set(handles.popupVel,'Enable','off');
    
    axes(handles.axesvelp)
    plot(Array.uMag(:,svvel)./nanmax(Array.uMag(:,svvel)),...
     Array.Depth_uMag(:,svvel)./Array.Bed(1,svvel),'ok')
    hold on     
    set(handles.axesvelp,'ydir','reverse');
    xlabel('u/umax','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2vel,'Enable','off');
    set(handles.r2vel,'String','');

    
elseif (strcmp(VelExtraM,'CON'))
    % Constant Extrapolation
    set(handles.popupVel,'Value',2);
    set(handles.popupVel,'Enable','off');

    axes(handles.axesvelp)
    plot(Array.uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
     Array.Depth_uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./Array.Bed(1,svvel),'ok')
    hold on     
    set(handles.axesvelp,'ydir','reverse');

    plot(Array.uMag(1:Array.indexUsurf(1,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
     Array.Depth_uMag(1:Array.indexUsurf(1,svvel),svvel)./Array.Bed(1,svvel),'oy','MarkerSize',9,'MarkerFaceColor','y')

    plot(Array.uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
     Array.Depth_uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./Array.Bed(1,svvel),'ob','MarkerSize',9,'MarkerFaceColor','b')

    xlabel('u/umax','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone','Surface Extrapolation','Bottom Extrapolation');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2vel,'Enable','off');
    set(handles.r2vel,'String','');
       
    
elseif (strcmp(VelExtraM,'LVI'))
    %3PT Slope
    set(handles.popupVel,'Value',2);
    set(handles.popupVel,'Enable','off');
    
    axes(handles.axesvelp)
    
     plot(Array.uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
         Array.Depth_uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./Array.Bed(1,svvel),'ok')
     hold on     
     set(handles.axesvelp,'ydir','reverse');
     
     plot(VelExtra.predsurf(:,svvel)./nanmax(Array.uMag(:,svvel)./100),...
         (Array.Bed(1,svvel)-C_vel.zpredsurfvel(:,svvel))./Array.Bed(1,svvel),'-g')
     
     plot(VelExtra.predbottom(:,svvel)./nanmax(Array.uMag(:,svvel)./100),...
         (Array.Bed(1,svvel)-C_vel.zpredbottomvel(:,svvel))./Array.Bed(1,svvel),'-g')

     plot(Array.uMag(1:Array.indexUsurf(1,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
         Array.Depth_uMag(1:Array.indexUsurf(1,svvel),svvel)./Array.Bed(1,svvel),'oy','MarkerSize',9,'MarkerFaceColor','y')

     plot(Array.uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
         Array.Depth_uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./Array.Bed(1,svvel),'ob','MarkerSize',9,'MarkerFaceColor','b')

    xlabel('u/umax','Fontsize',12);ylabel('z/H','Fontsize',12)    
    leg=legend('Measured Zone','Fitting Surf','Fitting Bot','Surface Extrapolation','Bottom Extrapolation');
    set(leg,'FontSize',9,'Location','best');
    ylim([0 1]);    xlim([0 1]);
    hold off
    set(handles.r2vel,'Enable','off');
    set(handles.r2vel,'String',''); 
            
elseif (strcmp(VelExtraM,'LPVI'))
    %Law of the Wall
    axes(handles.axesvelp)
    set(handles.popupVel,'Enable','on');
    sel=get(handles.popupVel,'Value');

    switch sel
        case 1%dimensional graph
             semilogx(Array.Bed(1,svvel)-Array.Depth_uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel),...
                 Array.uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)/100,'ok')
             hold on     

             if USys==0
                semilogx(VelExtra.zpred(:,svvel),VelExtra.pred(:,svvel),'-g')
             elseif USys==1
                semilogx(VelExtra.zpred(:,svvel)/0.3048,VelExtra.pred(:,svvel)/0.3048,'-g')
             end
             semilogx(Array.Bed(1,svvel)-Array.Depth_uMag(1:Array.indexUsurf(1,svvel),svvel),...
                 Array.uMag(1:Array.indexUsurf(1,svvel),svvel)/100,'oy','MarkerSize',9,'MarkerFaceColor','y')
             semilogx(Array.Bed(1,svvel)-Array.Depth_uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel),...
                 Array.uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)/100,'ob','MarkerSize',9,'MarkerFaceColor','b')

            xlabel('ln(z)','Fontsize',12);ylabel('u','Fontsize',12)        
            leg=legend('Measured Zone','Fitting','Surface Extrapolation','Bottom Extrapolation');
            set(leg,'FontSize',9,'Location','best');
            hold off
            
            %show r2 
            set(handles.r2vel,'Enable','on');
            if VelExtra.rcuad(1,svvel)<0.7
                set(handles.r2vel,'String',num2str(round(VelExtra.rcuad(1,svvel)*100)/100));
                set(handles.r2vel,'BackgroundColor','red');
            else
                set(handles.r2vel,'String',num2str(round(VelExtra.rcuad(1,svvel)*100)/100));
                set(handles.r2vel,'BackgroundColor','white');
            end
    
        case 2%adimensional graph

             plot(Array.uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
                 Array.Depth_uMag(Array.indexUmeas(1,svvel):Array.indexUmeas(2,svvel),svvel)./Array.Bed(1,svvel),'ok')
             hold on     
             set(handles.axesvelp,'ydir','reverse');
             plot(VelExtra.pred(:,svvel)./nanmax(Array.uMag(:,svvel)./100),...
                 (Array.Bed(1,svvel)-VelExtra.zpred(:,svvel))./Array.Bed(1,svvel),'-g')

             plot(Array.uMag(1:Array.indexUsurf(1,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
                 Array.Depth_uMag(1:Array.indexUsurf(1,svvel),svvel)./Array.Bed(1,svvel),'oy','MarkerSize',9,'MarkerFaceColor','y')

             plot(Array.uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./nanmax(Array.uMag(:,svvel)),...
                 Array.Depth_uMag(Array.indexUbottom(1,svvel):Array.indexUbottom(2,svvel),svvel)./Array.Bed(1,svvel),'ob','MarkerSize',9,'MarkerFaceColor','b')

            xlabel('u/umax','Fontsize',12);ylabel('z/H','Fontsize',12)    
            leg=legend('Measured Zone','Fitting','Surface Extrapolation','Bottom Extrapolation');
            set(leg,'FontSize',9,'Location','best');
            ylim([0 1]);    xlim([0 1]);
            hold off
            set(handles.r2vel,'Enable','off');
            set(handles.r2vel,'String',''); 
    
    end
    


end
    title('Velocity Profile','Fontsize',14)
  
