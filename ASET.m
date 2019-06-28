function varargout = ASET(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- ACOUSTIC SEDIMENT ESTIMATION TOOLBOX---
%
% ASET is a Matlab®-based toolbox developed to obtain, visualize, and 
% validate calibrations between the echo intensity level from static ADCP 
% measurements and suspended sediment concentration values measured 
% simultaneously in different points of surveyed verticals with traditional 
% techniques. In addition, once the calibration is obtained, ASET 
% transforms the acoustic intensity signal into sediment concentration for 
% each cell measured by the ADCP, obtaining the same spatial resolution as 
% the velocity field in the cross-section or temporal series (for static 
% measurements). Extrapolation methods are also included to estimate 
% velocity and concentration values in areas not measured by the ADCP, 
% (i.e., near bottom, surface, and banks), with visualization modules for 
% the user to evaluate each method. Finally, an integration module makes 
% possible to compute the total suspended sediment concentration in a 
% cross-section.

% CITATION:

%--------------------------------------------------------------------------
% Universidad Nacional del Litoral, Santa Fe, Argentina
% CONICET, Argentina
% U.S. Geological Survey
%
% Code contributed by Dominguez Ruben, L., Szupiany, R., Latosinski, F., 
% Lopez Weibel, C., Wood, M., Boldt, J.
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

try

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
catch err
    if isdeployed
        errLogFileName = fullfile(pwd,...
            ['errorLog' datestr(now,'yyyymmddHHMMSS') '.txt']);
        msgbox({['An unexpected error occurred. Error code: ' err.identifier];...
            ['Error details are being written to the following file: '];...
            errLogFileName},...
            'ASEt Status: Unexpected Error',...
            'error');
        fid = fopen(errLogFileName,'W');
        fwrite(fid,err.getReport('extended','hyperlinks','off'));
        fclose(fid);
        rethrow(err);        
    else
        close force
        msgbox(['An unexpected error occurred. Error code: ' err.identifier],...
            'ASET Status: Unexpected Error',...
            'error');
        rethrow(err);
    end
end

%add other folders
addpath 'tools'

warning off
% End initialization code - DO NOT EDIT


% --- Executes just before ASET is made visible.
function ASET_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ASET
handles.output = hObject;

% Read variables init
handles.readfilepd0=0; %pd0 default
handles.readadcpdata=0;%adcpdata default 
handles.VelExtraM='NONE';% Velocity extrapolation default
handles.CssExtraM='NONE';% Sediment extrapolation default
handles.units=0;% 0 metric (default); 1 english units
%Storage data
setappdata(0, 'units', handles.units)

guidata(hObject, handles);

% Initialize the GUI parameters:
% ------------------------------
guiparams.ASET_version = 'v1.5';

set(handles.figure1,'Name',['Acoustic Sediment Estimation Toolbox (ASET) ' guiparams.ASET_version], ...
    'DockControls','off')

% Draw the aset Background
% -----------------
pos = get(handles.figure1,'position');
axes(handles.ASETBackground);

 X = imread('ASET_Background.png');
 imdisp(X,'size',[pos(4) pos(3)]) % Avoids problems with users not having Image Processing TB

% Store the application ASET data:
% ---------------------------
setappdata(handles.figure1,'guiparams',guiparams)

% Initialize the GUI:
% -------------------
%initGUI(handles)
set_enable(handles,'init')
results_enable(handles,'init')
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ASET_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% Empty

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
%Empty

% --- Executes when user attempts to close_function figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

% Hint: delete(hObject) closes the figure
close_button = questdlg(...
    'You are about to exit ASET. Any unsaved work will be lost. Are you sure?',...
    'Exit ASET?','No');

switch close_button
    case 'Yes'
       close all force
       close all hidden
    otherwise
        return
end

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%
%Menu Toolbars
%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function file_function_Callback(hObject, ~, handles)
%empty

% --------------------------------------------------------------------
function new_function_Callback(hObject, eventdata, handles)
%Delete variables for new_function project
set_enable(handles,'init')
results_enable(handles,'init')
handles.VelExtraM='NONE';% Velocity extrapolation default
handles.CssExtraM='NONE';% Sediment extrapolation default
guidata(hObject,handles)

% Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%----------- ' datestr(now) ' ------------%'];...
        'New Project'};
    statusLogging(handles.LogWindow, log_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Open
% --------------------------------------------------------------------
function open_function_Callback(hObject, eventdata, handles)
% Empty
 

% --------------------------------------------------------------------
function singleselect_Callback(hObject, eventdata, handles)
%enable
set_enable(handles,'init')

persistent lastPath 
% If this is the first time running the function this session,
% Initialize lastPath to 0
if isempty(lastPath) 
    lastPath = 0;
end

%Start code
if lastPath == 0
[file,path] = uigetfile({'*ASC.TXT;*.mat',...
    'ASET Files (*.TXT,*.mat)';'*.*',  'All Files (*.*)'},'Select Input File');
else %remember the lastpath
[file,path] = uigetfile({'*ASC.TXT;*.mat',...
    'ASET Files (*.TXT,*.mat)';'*.*',  'All Files (*.*)'},'Select Input File',lastPath);
end

if file==0
else
    [guiparams]=ASET_ReadInputFile(path,file,lastPath);

    guiparams.R.multi=0; 
    guiparams.length=1;

    % Use the path to the last selected file
    % If 'uigetfile' is called, but no item is selected, 'lastPath' is not overwritten with 0
    if guiparams.PathName ~= 0
        lastPath = guiparams.PathName;
    end

    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
                '';...
                ['%--- ' datestr(now) ' ---%'];...
                'Input File Loaded';[cell2mat({guiparams.FileName})];
                'PD0 File Loaded';[cell2mat({guiparams.FileNamePD0})]};
                statusLogging(handles.LogWindow, log_text)
                
    %ADCP data
    [ADCPData]=ASET_ADCPData(guiparams);

    % Re-store the Application Data:
    % ------------------------------
    setappdata(0,'guiparams',guiparams)

    % enable boxes
     set_enable(handles,'loadfiles')
end


% --------------------------------------------------------------------
function multiselect_Callback(hObject, eventdata, handles)
set_enable(handles,'init')

handles.multiselect=ASET_MultiSelect;%MultiSelect Files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export Files
% --------------------------------------------------------------------
function export_function_Callback(hObject, eventdata, handles)
%Empty

% --------------------------------------------------------------------
function tecplotexport_Callback(hObject, eventdata, handles)

% Get the Application data:
% -------------------------
guiparams = getappdata(0,'guiparams');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File Variables exports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[file,path] = uiputfile({'*.dat','Save Tecplot Files (*.dat)'},...
    'Save Tecplot Files');

if ischar(file)
    waitmessage = ['Exporting Tecplot File...' file];
    hwait = waitbar(0,waitmessage,'Name','ASET');
    
    guiparams.savefiletecplot=fullfile(path,file);

    % Re-store the Application data:
    % ------------------------------
    setappdata(0,'guiparams',guiparams)
        
     %Create the Tecplot variables file
     ASET_CreateTecplotFile_withExtrapo(guiparams,fullfile(path,file),handles);

     waitbar(1/2,hwait)
     %Create the Bathymetry Tecplot file
     ASET_CreateTecplotFileBat(guiparams,fullfile(path,file),handles);
        
     waitbar(1,hwait)
     delete(hwait)
     
    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'Complete Export for Tecplot file';fullfile(path,file)};
        statusLogging(handles.LogWindow, log_text)
    
end


% --------------------------------------------------------------------
function excelexport_Callback(hObject, eventdata, handles)
%Export EXCEL File
ASET_ExportXlsFile(handles);%  


% --------------------------------------------------------------------
function matfile_Callback(hObject, eventdata, handles)
%Export mat file
guiparams = getappdata(0,'guiparams');

[file,path] = uiputfile('*.mat','Save file');
if file==0
else
    hwait = waitbar(0,'Exporting .mat File...','Name','ASET');
    save([path file], 'guiparams');
    waitbar(1,hwait)
    delete(hwait)
    
    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'Complete Export .mat file' ;[path file]};
        statusLogging(handles.LogWindow, log_text)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphics
% --------------------------------------------------------------------
function graphics_function_Callback(hObject, eventdata, handles)
%empty


% --------------------------------------------------------------------
function plotcrosssection_Callback(hObject, eventdata, handles)
ControlVar(handles);

%Variables plot cross section
handles.getplotcross = ASET_ControlPlot;
uiresume(gcbf);

% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'Plot Cross Section Complete'};
statusLogging(handles.LogWindow, log_text)


% --------------------------------------------------------------------
function plotverticalprofile_Callback(hObject, eventdata, handles)
ControlVar(handles);
%Variables plot cross section
handles.getplotcross = ASET_VerticalProfilePlot;
uiresume(gcbf);

% Push messages to Log Window:
% ----------------------------
log_text = {...
'';...
['%--- ' datestr(now) ' ---%'];...
'Plot Vertical Profilers Complete'};
statusLogging(handles.LogWindow, log_text)
    
    
 % --------------------------------------------------------------------
function close_function_Callback(hObject, eventdata, handles)
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tools
% --------------------------------------------------------------------
function tool_function_Callback(hObject, eventdata, handles)
% empty


% --------------------------------------------------------------------
function plottest_Callback(hObject, eventdata, handles)
ControlVar(handles);

handles.getplottest = ASET_PlotTest;
uiresume(gcbf);


% --------------------------------------------------------------------
function calibrationmethodp_Callback(hObject, eventdata, handles)
handles.getcalibration = ASET_Calibration;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help
% --------------------------------------------------------------------
function help_function_Callback(hObject, eventdata, handles)
%empty


% --------------------------------------------------------------------
function usersguide_Callback(hObject, eventdata, handles)
try
web('https://asetoolbox.blogspot.com/p/blog-page_97.html')%CHANGE
catch err %#ok<NASGU>
	if isdeployed
        errLogFileName = fullfile(pwd,...
            ['errorLog' datestr(now,'yyyymmddHHMMSS') '.txt']);
        msgbox({['An unexpected error occurred. Error code: ' err.identifier];...
            ['Error details are being written to the following file: '];...
            errLogFileName},...
            'mstat Status: Unexpected Error',...
            'error');
        fid = fopen(errLogFileName,'W');
        fwrite(fid,err.getReport('extended','hyperlinks','off'));
        fclose(fid);
        rethrow(err)
    else
        msgbox(['An unexpected error occurred. Error code: ' err.identifier],...
            'mstat Status: Unexpected Error',...
            'error');
        rethrow(err);
    end
end


% --------------------------------------------------------------------
function updates_Callback(hObject, eventdata, handles)
try
web('https://asetoolbox.blogspot.com/p/blog-page_3.html')%CHANGE
catch err %#ok<NASGU>
	if isdeployed
        errLogFileName = fullfile(pwd,...
            ['errorLog' datestr(now,'yyyymmddHHMMSS') '.txt']);
        msgbox({['An unexpected error occurred. Error code: ' err.identifier];...
            ['Error details are being written to the following file: '];...
            errLogFileName},...
            'mstat Status: Unexpected Error',...
            'error');
        fid = fopen(errLogFileName,'W');
        fwrite(fid,err.getReport('extended','hyperlinks','off'));
        fclose(fid);
        rethrow(err)
    else
        msgbox(['An unexpected error occurred. Error code: ' err.identifier],...
            'mstat Status: Unexpected Error',...
            'error');
        rethrow(err);
    end
end


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
try
web('https://asetoolbox.blogspot.com/')%CHANGE
catch err %#ok<NASGU>
	if isdeployed
        errLogFileName = fullfile(pwd,...
            ['errorLog' datestr(now,'yyyymmddHHMMSS') '.txt']);
        msgbox({['An unexpected error occurred. Error code: ' err.identifier];...
            ['Error details are being written to the following file: '];...
            errLogFileName},...
            'mstat Status: Unexpected Error',...
            'error');
        fid = fopen(errLogFileName,'W');
        fwrite(fid,err.getReport('extended','hyperlinks','off'));
        fclose(fid);
        rethrow(err)
    else
        msgbox(['An unexpected error occurred. Error code: ' err.identifier],...
            'mstat Status: Unexpected Error',...
            'error');
        rethrow(err);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial Panel 

% --- Executes on button press in adcpdataoption.
function adcpdataoption_Callback(hObject, eventdata, handles)

ADCPdata=getappdata(0, 'ADCPdata');
ASET_ADCPData(ADCPdata)


% --- Executes on button press in fileoptions.
function fileoptions_Callback(hObject, eventdata, handles)
handles.fileoptions=ASET_FileOptions;


%%%%%%%%%%%%%%%%%%%
% Input Parameters 
%%%%%%%%%%%%%%%%%%%
function finediameter_Callback(hObject, eventdata, handles)
% Empty

% --- Executes during object creation, after setting all properties.
function finediameter_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function coursediameter_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function coursediameter_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fineconcentration_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function fineconcentration_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%
%Extrapolations Methods
%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in VelExtrapolation.
function VelExtrapolation_Callback(hObject, eventdata, handles)
sel=get(handles.VelExtrapolation,'Value');

switch sel
    case 1
        handles.VelExtraM='NONE';
        guidata(hObject,handles)
    case 2
        handles.VelExtraM='CON';
        guidata(hObject,handles)
    case 3
        handles.VelExtraM='LVI';%3Pt. Slope
        guidata(hObject,handles)
    case 4
        handles.VelExtraM='LPVI';%Law of the Wall
        guidata(hObject,handles)
end
guidata(hObject,handles)

 
% --- Executes during object creation, after setting all properties.
function VelExtrapolation_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CssExtrapolation.
function CssExtrapolation_Callback(hObject, eventdata, handles)

sel=get(handles.CssExtrapolation,'Value');

switch sel
    case 1
        handles.CssExtraM='NONE';
        guidata(hObject,handles)
    case 2
        handles.CssExtraM='CON';%Constant
        guidata(hObject,handles)
    case 3
        handles.CssExtraM='LSI';%3Pt. Slope
        guidata(hObject,handles)
    case 4
        handles.CssExtraM='RSI';%Rouse 
        guidata(hObject,handles)
end
 guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CssExtrapolation_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lateralextraoption.
function lateralextraoption_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function lateralextraoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lateralextraoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%

function surface_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of surface as text
%        str2double(get(hObject,'String')) returns contents of surface as a double


% --- Executes during object creation, after setting all properties.
function surface_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function measurement_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function measurement_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bottom_Callback(hObject, eventdata, handles)
% empty

% --- Executes during object creation, after setting all properties.
function bottom_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function total_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function total_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function gw_Callback(hObject, eventdata, handles)
% empty

% --- Executes during object creation, after setting all properties.
function gw_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function totaldischarge_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function totaldischarge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totaldischarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Calculate
function calculate_Callback(hObject, eventdata, handles)
%Calculate transport
set_enable(handles,'results')

%Control input var
Control=ControlVar(handles);

if Control==3; %all variables included
    % Get the Application data:
    % -------------------------
    guiparams = getappdata(0,'guiparams');
    VelExtraM   = handles.VelExtraM;%Interpolation of velocity
    CssExtraM  = handles.CssExtraM;%Interpolation of sedimento
    Inst= guiparams.Inst;%PD0 data
    Sensor= guiparams.Sensor;%PD0 data
    Cfg= guiparams.Cfg;%PD0 data
    LatExtra=get(handles.lateralextraoption,'Value');%Lateral extrapolation methods
    ADCPdata=getappdata(0, 'ADCPdata');

    % Re-store the Application Data:
    % ------------------------------
    setappdata(0,'readformatfile',guiparams.readformatfile)

    V   = guiparams.V;
    A   = guiparams.A;
    Edge = guiparams.Edge;
    X   = guiparams.xutm;
    Y   = guiparams.yutm;        
    Filename = guiparams.FileName;

    %Read data input fine material diameter [mm]
    finediameter=get(handles.finediameter,'String');
    S.Df=str2num(finediameter)/1000;

    %Read data input course material diameter [mm]
    coursediameter=get(handles.coursediameter,'String');
    S.Dc=str2num(coursediameter)/1000;

    %Read data input fine material concentration [mg/l]
    fineconcentration=get(handles.fineconcentration,'String');
    S.Csf=str2num(fineconcentration);

    %Course material concentration [mg/l]
    S.Css=10;%Default value
    
    waitmessage=['Processing data...' Filename];
    handles.hwait = waitbar(0,waitmessage,'Name','ASET',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
    guidata(hObject,handles) 
    setappdata(handles.hwait,'canceling',0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %Static Parameters
    M=1;% calculate of Gss
    Mont=0;%Without MonteCarlo
    R=[];%Without R variable

    [St]=ASET_ParamStatic(M,V);

    [P]=ASET_Parameters(M,V,S,R,Inst,Cfg,Sensor,St,ADCPdata,Mont,1);

    waitbar(0.1,handles.hwait);

    %
    [Cut_vel]=ASET_ArrayResizeVel(V,St);

    [Cut_back]=ASET_ArrayResizeBack(V,St);
    waitbar(0.4,handles.hwait);

    %
    [Css,~,~,P]=ASET_Concentration(M,V,S,Cut_back,R,St,P);
    waitbar(0.5,handles.hwait);

    %    
    if (strcmp(VelExtraM,'NONE')) & (strcmp(CssExtraM,'NONE')) | (strcmp(VelExtraM,'NONE')) %without velocity extrapolation
        VelExtra=[];
        C_vel=[];
        if (strcmp(VelExtraM,'NONE')) & (strcmp(CssExtraM,'NONE'))
            C_back=[];
            CssExtra=[];
        elseif (strcmp(VelExtraM,'NONE'))
            [C_back]=ASET_ArrayResizeBackExtrap(V,St,Cut_back);
            [CssExtra]=ASET_ExtrapCss(V,S,Css,Cut_back,St,CssExtraM,P,C_back);
        end
        [Dis]=ASET_DischargeLiquid(V,VelExtra,C_vel,St,A,guiparams.readformatfile);
        [Meas]=ASET_TransportMeas(Css,V,S,St,Dis);
        [NoMeas]=ASET_TransportNoMeas(V,S,CssExtra,VelExtra,C_back,C_vel,St,Dis,Cut_vel,Cut_back);
    elseif (strcmp(CssExtraM,'NONE'))%without concetration extrapolation
        [C_vel]=ASET_ArrayResizeVelExtrap(V,St,Cut_vel);
        [VelExtra]=ASET_ExtrapV(V,Cut_vel,St,VelExtraM,C_vel);
        C_back=[];
        CssExtra=[];
        [Dis]=ASET_DischargeLiquid(V,VelExtra,C_vel,St,A,guiparams.readformatfile);
        [Meas]=ASET_TransportMeas(Css,V,S,St,Dis);
        [NoMeas]=ASET_TransportNoMeas(V,S,CssExtra,VelExtra,C_back,C_vel,St,Dis,Cut_vel,Cut_back);    
    else %with extrapolations both variables
        [C_vel]=ASET_ArrayResizeVelExtrap(V,St,Cut_vel);
        [C_back]=ASET_ArrayResizeBackExtrap(V,St,Cut_back);
        [VelExtra]=ASET_ExtrapV(V,Cut_vel,St,VelExtraM,C_vel);
        [CssExtra]=ASET_ExtrapCss(V,S,Css,Cut_back,St,CssExtraM,P,C_back);
        [Dis]=ASET_DischargeLiquid(V,VelExtra,C_vel,St,A,guiparams.readformatfile);
        [Meas]=ASET_TransportMeas(Css,V,S,St,Dis);
        [NoMeas]=ASET_TransportNoMeas(V,S,CssExtra,VelExtra,C_back,C_vel,St,Dis,Cut_vel,Cut_back);
    end

    %
    waitbar(0.6,handles.hwait);

    %concatenate all data
    [Array]=ASET_Concatenate(V,Css,Cut_vel,Cut_back,VelExtra,CssExtra,St,Meas,NoMeas,C_vel,C_back);

    %Lateral Extrapolation 
    if LatExtra==1%without lateral extrapolation
        QEdge.Gss.left=0;
        QEdge.Gss.right=0;
        QEdge.Gw.left=0;
        QEdge.Gw.right=0;
        QEdge.Q.right=0;
        QEdge.Q.left=0;
    else
        if LatExtra==2%triangular method
            Method.Edge=0;%number of method used on ASET_LateralExtraDis
        elseif LatExtra==3%rectangular method
            Method.Edge=1;%number of method used on ASET_LateralExtraDis
        end

        [QEdge]=ASET_LateralExtraDis(Method,Array,V,S,Edge);
    end

    %Final the calculate

      %enable nomeasured zones
         set(handles.bottom,'Enable','on');
         set(handles.surface,'Enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coarse material transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NoMeas.Surf.GssT=nansum(nansum(NoMeas.Surf.coarse));%Surface zone transport
        NoMeas.Bottom.GssT=nansum(nansum(NoMeas.Bottom.coarse));%Bottom zone trasport
        
        set(handles.surface,'String',num2str(round(NoMeas.Surf.GssT,2)));
        set(handles.bottom,'String',num2str(round(NoMeas.Bottom.GssT,2)));
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fine material transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NoMeas.Surf.GwT=nansum(nansum(NoMeas.Surf.fine));%Surface zone transport
        NoMeas.Bottom.GwT=nansum(nansum(NoMeas.Bottom.fine));%Bottom zone trasport

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coarse material transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Meas.GssT=nansum(nansum(Meas.coarse));%Measurement zone trasport
        Array.GssT=QEdge.Gss.left+QEdge.Gss.right+NoMeas.Surf.GssT+...
            Meas.GssT+NoMeas.Bottom.GssT;%Total transport coarse material

        %write results coarse material
        set(handles.measurement,'String',num2str(round(Meas.GssT,2)));
        set(handles.total,'String',num2str(round(Array.GssT,2)));
        set(handles.totaldischarge,'String',num2str(round(Dis.Total+QEdge.Q.right+QEdge.Q.left,2)));         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fine material transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Meas.GwT=nansum(nansum(Meas.fine));%Measurement zone trasport
    Array.GwT=NoMeas.Surf.GwT+Meas.GwT+NoMeas.Bottom.GwT+...
        QEdge.Gw.left+QEdge.Gw.right;%fine material
       
    %Write results fine material
     set(handles.gw,'String',num2str(round(Array.GwT,2)));%Total fine material trasport
        

    % Push messages to Log Window:
    % ----------------------------
    log_text = {
                ['%--- ' datestr(now) ' ---%'];...
                sprintf(' File analyzed: %s', guiparams.FileName);...
                sprintf(' Navigation Reference: %s',cell2mat(guiparams.V.navref));...
                sprintf(' Distance cross section [m]: %0.1f',V.mcsDist(1,end));...
                sprintf(' Average depth [m]: %0.1f',nanmean(Array.Bed));...
                sprintf(' Average velocity [cm/s]: %0.2f',nanmean(nanmean(Array.uMag)));...
                sprintf(' Q surface [m3/s]: %0.1f',Dis.surfT);...
                sprintf(' Q bottom [m3/s]: %0.1f',Dis.bottomT);...
                sprintf(' Q left [m3/s]: %0.2f',QEdge.Q.left);...
                sprintf(' Q right [m3/s]: %0.2f',QEdge.Q.right);...
                sprintf(' Average Ms2 [mg/l]: %0.1f',nanmean(nanmean(Css))*1000);...
                sprintf(' Qss left [kg/s]: %0.2f',QEdge.Gss.left);...
                sprintf(' Qss right [kg/s]: %0.2f',QEdge.Gss.right)};
                statusLogging(handles.LogWindow, log_text)

                
    waitbar(1,handles.hwait)
    delete(handles.hwait)

    %Storage data
    guiparams.S=S;
    guiparams.V=V;
    guiparams.Dis=Dis;
    guiparams.VelExtraM=VelExtraM;
    guiparams.CssExtraM=CssExtraM;
    guiparams.LatExtra=LatExtra;
    guiparams.St=St;
    guiparams.P=P;
    guiparams.C_vel=C_vel;
    guiparams.C_back=C_back;
    guiparams.Css=Css;
    guiparams.Meas=Meas;
    guiparams.Array=Array;
    guiparams.Pos.xutm=X;
    guiparams.Pos.yutm=Y;
    guiparams.QEdge=QEdge;
    guiparams.NoMeas=NoMeas;
    guiparams.VelExtra=VelExtra;
    guiparams.CssExtra=CssExtra;


    %Storage data
    setappdata(0, 'guiparams', guiparams)

    set(handles.measurement,'Enable','on');
    set(handles.total,'Enable','on');
    set(handles.gw,'Enable','on');
    set(handles.totaldischarge,'Enable','on');
    results_enable(handles,'loadfiles')


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Graphics results profiles
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    cof.horizontalsmooth=0;

    cof.verticalsmooth=0;

    if (strcmp(CssExtraM,'NONE')) 
        Array.Css=smooth2a(Css,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_Css=smooth2a(V.mcsDepth,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(V.mcsDist,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Bed=V.mcsBed;
    else
        Array.Css=smooth2a(Array.Css,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Depth_Css=smooth2a(Array.Depth_Css,cof.horizontalsmooth,cof.verticalsmooth);
        Array.Dist=smooth2a(Array.Dist,cof.horizontalsmooth,cof.verticalsmooth);
    end

    axes(handles.crosssectionaxes)

    pcolor(Array.Dist,Array.Depth_Css,Array.Css*1000);
    shading interp

    title('ASET: Mean Cross Section Contour Ms2 [mg/l]','Fontsize',10)
    xlabel('Distance across section [m]','Fontsize',8);ylabel('Depth [m]','Fontsize',8);
    fig.Label.String = 'Ms2 Sand Concentration [mg/l]';


    grid off
    hold on
    plot(Array.Dist(1,1:end),Array.Bed,'-k','Linewidth',2)
    hold off
    set(gca,'ydir','reverse');
    ylim([0 max(max(Array.Bed))]);
    fig.xlim=get(gca,'xlim');
    fig.ylim=get(gca,'ylim');
    fig.zmin=1.35*floor(nanmin(nanmin(Array.Css*1000)));
    fig.zmax=0.25*ceil(nanmax(nanmax(Array.Css*1000)));
    caxis([fig.zmin fig.zmax]);
    fig.colormap=colormap('jet');
    fig.colorbar=colorbar;
    fig.FontSize=16;

    %%
    %Vertical profiles

    %Vertical profile Concentration
    axes(handles.verticalprofileaxes)
    plot(Array.Css_ave,Array.Depth_Css_ave,'o','MarkerEdgeColor',[0.5 0.5 0.5],'Markersize',1);
    hold on
    set(gca,'ydir','reverse')
    ylim([0 1]);xlim([0 1])   
    title('ASET: Ms2 Vertical Profile','Fontsize',8)
    ylabel('Relative Depth','Fontsize',8);xlabel('Ms2/Ms2max','Fontsize',8)
    hold off

    %Vertical profile Velocity
    axes(handles.verticalprofileaxes1)
    plot(Array.uMag_ave,Array.Depth_uMag_ave,'o','MarkerEdgeColor',[0.5 0.5 0.5],'Markersize',1);
    hold on
    set(gca,'ydir','reverse')
    ylim([0 1]);xlim([0 1])
    title('ASET: Velocity Vertical Profile','Fontsize',8)
    ylabel('Relative Depth','Fontsize',8);xlabel('u/umax','Fontsize',8)
    hold off

end
    

% --- Executes on selection change in LogWindow.
function LogWindow_Callback(hObject, eventdata, handles)
%Empty


% --- Executes during object creation, after setting all properties.
function LogWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LogWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_enable(handles,enable_state)

switch enable_state
    case 'init'
        set(handles.CssExtrapolation,'Value',1);
        set(handles.VelExtrapolation,'Value',1);
        set(handles.lateralextraoption,'Value',1);
        set(handles.surface,'String','','Enable','off');
        set(handles.measurement,'String','','Enable','off');
        set(handles.bottom,'String','','Enable','off');
        set(handles.total,'String','','Enable','off');
        set(handles.gw,'String','','Enable','off');
        set(handles.VelExtrapolation,'Enable','off');
        set(handles.CssExtrapolation,'Enable','off');
        set(handles.lateralextraoption,'Enable','off');
        set(handles.finediameter,'String','','Enable','off');
        set(handles.coursediameter,'String','','Enable','off');
        set(handles.fineconcentration,'String','','Enable','off');
        set(handles.calculate,'Enable','off');
        set(handles.adcpdataoption,'Enable','off');
        set(handles.multiselect,'Enable','on');
        set(handles.totaldischarge,'String','','Enable','off');
        set(handles.plottest,'Enable','off');
        cla(handles.crosssectionaxes,'reset')
        axes(handles.crosssectionaxes)
        grid on
        cla(handles.verticalprofileaxes,'reset')
        axes(handles.verticalprofileaxes)
        grid on
        cla(handles.verticalprofileaxes1,'reset')
        axes(handles.verticalprofileaxes1)
        grid on
    case 'loadfiles'
        set(handles.surface,'Enable','off');
        set(handles.measurement,'Enable','off');
        set(handles.bottom,'Enable','off');
        set(handles.total,'Enable','off');
        set(handles.gw,'Enable','off');       
        set(handles.VelExtrapolation,'Enable','on');
        set(handles.CssExtrapolation,'Enable','on');
        set(handles.finediameter,'Enable','on');
        set(handles.coursediameter,'Enable','on');
        set(handles.fineconcentration,'Enable','on');
        set(handles.lateralextraoption,'Enable','on');
        set(handles.calculate,'Enable','on');    
        set(handles.adcpdataoption,'Enable','on');
        set(handles.multiselect,'Enable','on');
        set(handles.totaldischarge,'Enable','off');
    case 'results'
        cla(handles.crosssectionaxes,'reset')
        cla(handles.verticalprofileaxes,'reset')
        cla(handles.verticalprofileaxes1,'reset')
        set(handles.surface,'String','','Enable','off');
        set(handles.measurement,'String','','Enable','off');
        set(handles.bottom,'String','','Enable','off');
        set(handles.total,'String','','Enable','off');
        set(handles.gw,'String','','Enable','off');
        set(handles.totaldischarge,'String','','Enable','off');
    otherwise
end


% 
function results_enable(handles,enable_state)
switch enable_state
    case 'init'
        set(handles.plottest,'Enable','off');
        set(handles.export_function,'Enable','off');
    case 'loadfiles'
        set(handles.export_function,'Enable','on');
        set(handles.plottest,'Enable','on');
    otherwise
end
            

function [Control]=ControlVar(handles)

if isempty(get(handles.finediameter,'String')) 
    wrong_data=warndlg('Please complete the Fine Sediment Diameter box',...
    'WARNING'); 
    Control1=0;
else
    Control1=1;
end
if isempty(get(handles.coursediameter,'String'))
    wrong_data=warndlg('Please complete the Coarse Sediment Diameter box',...
    'WARNING');
    Control2=0;
else
    Control2=1;
end
if isempty(get(handles.fineconcentration,'String'))
    wrong_data=warndlg('Please complete the Fine Sediment Concentration box',...
    'WARNING'); 
    Control3=0;
else
    Control3=1;
end

    Control=Control1+Control2+Control3;
