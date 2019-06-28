function varargout = ASET_MultiSelect(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function read and calculate more than one file. The methodology of
% calculate is the same of the initial panel.

%by Dominguez Ruben,L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_MultiSelect_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_MultiSelect_OutputFcn, ...
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


% --- Executes just before ASET_MultiSelect is made visible.
function ASET_MultiSelect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for ASET_MultiSelect
handles.output = hObject;
handles.excelimport=0;
handles.units=0;%default by SI system
handles.quitfile=0;
handles.first=1;
% Update handles structure
guidata(hObject, handles);

%initGUI
set_enable(handles,'init')

set(handles.figure1,'Name',['ASET: Multiselect '], ...
    'DockControls','off')


% UIWAIT makes ASET_MultiSelect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ASET_MultiSelect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%Toolbar menu
%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function filefunctions_Callback(hObject, eventdata, handles)
% Empty


% --------------------------------------------------------------------
function newmultiselect_Callback(hObject, eventdata, handles)
set_enable(handles,'init')
handles.first=1;
guidata(hObject, handles)
% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%----------- ' datestr(now) ' ------------%'];...
    'New Project'};
statusLogging(handles.LogWindow, log_text)

% --------------------------------------------------------------------
function toolbar_Callback(hObject, eventdata, handles)
% Empty


% --------------------------------------------------------------------
function exportfunction_Callback(hObject, eventdata, handles)
% Empty


% --------------------------------------------------------------------
function exporttecplot_Callback(hObject, eventdata, handles)
%this function export data in .dat format for Tecplot
    %Input data
    guiparams = getappdata(0,'guiparams');
for i=1:handles.numfile

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %File Variables exports
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    message=['Save Tecplot Files of ' guiparams{i}.FileName];
    [file,path] = uiputfile({'*.dat','Save Tecplot Files (*.dat)'},...
        message);


    if ischar(file)
        waitmessage = ['Exporting Tecplot File...' file];
        hwait = waitbar(0,waitmessage,'Name','ASET');

        guiparams{i}.savefiletecplot=fullfile(path,file);

         %Create the Tecplot variables file
         %with extrapolate zone
         Vartec=guiparams{i};
         ASET_CreateTecplotFile_withExtrapo(Vartec,fullfile(path,file),handles);
        

         waitbar(1/2,hwait)
         %Create the Bathymetry Tecplot file
         ASET_CreateTecplotFileBat(Vartec,fullfile(path,file),handles);
         
         % Push messages to Log Window:
        % ----------------------------
        log_text = {...
            '';...
            ['%--- ' datestr(now) ' ---%'];...
            'Complete Export .dat file';[path file]};
        statusLogging(handles.LogWindow, log_text)
        
         waitbar(1,hwait)
         delete(hwait)
    else
        break
    end

clear Array Vartec
end



% --------------------------------------------------------------------
function exportmatlab_Callback(hObject, eventdata, handles)
%this funtion exporta all data on .mat format
hwait = waitbar(0,'Exporting .mat File...','Name','ASET');

guiparams = getappdata(0,'guiparams');

[file,path] = uiputfile('*.mat','Save MATLAB File');%name and path

if file==0
    %cancel save
    delete(hwait)
else
    save([path file], 'guiparams');

    waitbar(1,hwait)
    delete(hwait)

    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'Complete Export .mat File' ;[path file]};
    statusLogging(handles.LogWindow, log_text)
end


% --------------------------------------------------------------------
function exportexcel_Callback(hObject, eventdata, handles)
%Export excel data
ASET_ExportXlsFile(handles);%  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual cross section or measurements
% --------------------------------------------------------------------

% --- Executes on button press in addfiles.
function addfiles_Callback(hObject, eventdata, handles)


persistent lastPath 
% If this is the first time running the function this session,
% Initialize lastPath to 0
if isempty(lastPath) 
    lastPath = 0;
end

%Start code
if lastPath == 0
[file,path] = uigetfile({'*ASC.TXT;*.mat',...
    'ASET Files (*.TXT,*.mat)';'*.*',  'All Files (*.*)'},'Select Input File','MultiSelect','on');
else %remember the lastpath
[file,path] = uigetfile({'*ASC.TXT;*.mat',...
    'ASET Files (*.TXT,*.mat)';'*.*',  'All Files (*.*)'},'Select Input File','MultiSelect','on',lastPath);
end

if handles.first==1
    %enable
    set_enable(handles,'init')
    if iscell(file)
        for i=1:length(file)
            [guiparams{i}]=ASET_ReadInputFile(path,file{i},lastPath);

            guiparams{i}.R.multi=1; 
            guiparams{i}.length=length(file);

            % Use the path to the last selected file
            % If 'uigetfile' is called, but no item is selected, 'lastPath' is not overwritten with 0
            if guiparams{i}.PathName ~= 0
                lastPath = guiparams{i}.PathName;
            end

            % Push messages to Log Window:
            % ----------------------------
            log_text = {...
                        '';...
                        ['%--- ' datestr(now) ' ---%'];...
                        'Input File Loaded';[cell2mat({guiparams{i}.FileName})];
                        'PD0 File Loaded';[cell2mat({guiparams{i}.FileNamePD0})]};
                        statusLogging(handles.LogWindow, log_text)
        end

        handles.numfile=length(file);

        %ADCP data
        [ADCPData]=ASET_ADCPData(guiparams{1});

        % Re-store the Application Data:
        % ------------------------------
        setappdata(0,'guiparams',guiparams)

        %Show data input files table
        handles.fileinput=cell(handles.numfile,2);
        handles.fileinput(:,1)=file;
        for i=1:handles.numfile
            handles.fileinput(i,2)={guiparams{i}.FileNamePD0};
        end
        guidata(hObject,handles) 
        set(handles.inputfiles,'Data',handles.fileinput)

        %Show data table 2. Sedimentology table
        celltable1=cell(handles.numfile,6);
        celltable1(:,1:3)={''};
        celltable1(:,4:6)={'None'};
        set(handles.sedimentaTable,'Data',celltable1)

        %set enable
        set_enable(handles,'loadfiles')

        % enable boxes
        set_enable(handles,'loadfiles')
        handles.first=0;
        guidata(hObject,handles) 
    else
    end
elseif handles.first==0
    if ischar(file)
        handles.length=handles.numfile+size(file,1);
    elseif iscell(file)
        handles.length=handles.numfile+length(file);
    end
    y=1;
        for i=handles.numfile+1:handles.length
            guiparams = getappdata(0,'guiparams');
            
            if ischar(file)
                [guiparams{i}]=ASET_ReadInputFile(path,file,lastPath);
            elseif iscell(file)
                [guiparams{i}]=ASET_ReadInputFile(path,file{y},lastPath);
            end
            
            guiparams{i}.R.multi=1; 
            guiparams{i}.length=length(file);

            % Use the path to the last selected file
            % If 'uigetfile' is called, but no item is selected, 'lastPath' is not overwritten with 0
            if guiparams{i}.PathName ~= 0
                lastPath = guiparams{i}.PathName;
            end

            % Push messages to Log Window:
            % ----------------------------
            log_text = {...
                        '';...
                        ['%--- ' datestr(now) ' ---%'];...
                        'Input File Loaded';[cell2mat({guiparams{i}.FileName})];
                        'PD0 File Loaded';[cell2mat({guiparams{i}.FileNamePD0})]};
                        statusLogging(handles.LogWindow, log_text)

                        y=y+1;
        end

        %handles.numfile=handles.length;

        %ADCP data
        [ADCPData]=ASET_ADCPData(guiparams{1});

        % Re-store the Application Data:
        % ------------------------------
        setappdata(0,'guiparams',guiparams)
        
        TableDATA = get(handles.inputfiles,'Data');
        TableDATA(handles.numfile+1:handles.length,1) = {file};
         for i=handles.numfile+1:handles.length
            TableDATA(i,2) ={guiparams{i}.FileNamePD0};
         end
        set(handles.inputfiles,'Data',TableDATA)
        
        TableDATASED = get(handles.sedimentaTable,'Data');
        %Show data table 2. Sedimentology table
        celltable1=cell(handles.length,6);
        celltable1(:,1:3)={''};
        celltable1(:,4:6)={'None'};
        set(handles.sedimentaTable,'Data',celltable1)    

        % enable boxes
        set_enable(handles,'loadfiles')
        
        handles.numfile=handles.length;
        guidata(hObject,handles)
end



% --- Executes on button press in quitfiles.
function quitfiles_Callback(hObject, eventdata, handles)
if handles.quitfile==0
   warndlg('Please select the file that you can delete')
else
    %delete the data selected 
    TableDATA = get(handles.inputfiles,'Data');

    TableDATASED = get(handles.sedimentaTable,'Data');

    TableDATA(handles.quitfile,:)=[];

    set(handles.inputfiles,'Data',TableDATA)

    TableDATASED(handles.quitfile,:)=[];
    set(handles.sedimentaTable,'Data',TableDATASED)

    guiparams = getappdata(0,'guiparams');

    guiparams(handles.quitfile)=[];
    
    % Re-store the Application Data:
    % ------------------------------
    setappdata(0,'guiparams',guiparams)

    handles.numfile=handles.numfile-length(handles.quitfile);

    guidata(hObject,handles)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%Main Panel
%%%%%%%%%%%%%%%%%%
     
% --- Executes on button press in adcpdataoption.
function adcpdataoption_Callback(hObject, eventdata, handles)
%ADCP data
ADCPdata=getappdata(0, 'ADCPdata');
ASET_ADCPData(ADCPdata)


% --- Executes on button press in filesdata.
function filesdata_Callback(hObject, eventdata, handles)
handles.fileoptions=ASET_FileOptions;


% --- Executes on button press in importdata.
function importdata_Callback(hObject, eventdata, handles)

%edit table (empty)
set(handles.sedimentaTable, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) size(get(handles.sedimentaTable,'ColumnFormat'),2)]));
set(handles.results, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) size(get(handles.results,'Data'),2)]));
set(handles.SumT, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) 11]));
    
%Import Sedimentology data from Excel files
[filenamexls,pathnamexls] = uigetfile({'*.xlsx';'*.*'},'Select *.xlsx file');
 
if ischar(filenamexls)
           
    waitmessage=['Importing data...' filenamexls];
    handles.hwait = waitbar(0,waitmessage,'Name','ASET');
    infile=[pathnamexls filenamexls];  
    xls=xlsread(infile);
    waitbar(0.5,handles.hwait);

   if ischar(xls(1,1))
       xls(1,:)=[];
   else
   end
   
   if size(xls,1)==handles.numfile

    for t=1:size(xls)
    %Velocity method
        if xls(t,4)==1
           VelExtraM{t,1}='None';
        elseif xls(t,4)==2
           VelExtraM{t,1}='Constant';%Linear Velocity 
        elseif xls(t,4)==3
           VelExtraM{t,1}='3Pt Slope';%Linear Velocity 
        elseif xls(t,4)==4
           VelExtraM{t,1}='Law of the Wall';%Log Profile velocity interpolation
        end
       waitbar(0.6,handles.hwait);
    %Concetration method
        if xls(t,5)==1
            CssExtraM{t,1}='None';
        elseif xls(t,5)==2
            CssExtraM{t,1}='Constant';%Linear Method
        elseif xls(t,5)==3
            CssExtraM{t,1}='3Pt Slope';%Linear Method
        elseif xls(t,5)==4
            CssExtraM{t,1}='Rouse';%Rouse Method
        end
       waitbar(0.7,handles.hwait);
    %Lateral extrapolation
        if xls(t,6)==1
            LatExtra{t,1}='None';
        elseif xls(t,6)==2
            LatExtra{t,1}='Triangular';
        elseif xls(t,6)==3
            LatExtra{t,1}='Rectangular';
        end
    waitbar(0.8,handles.hwait);

    end

    set(handles.sedimentaTable, 'data',[num2cell(xls(:,1:3)) VelExtraM CssExtraM LatExtra]);

    handles.excelimport=1;
    guidata(hObject,handles)

    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'Import Excel File Complete'};
    statusLogging(handles.LogWindow, log_text)
        waitbar(1,handles.hwait);
        delete(handles.hwait)
        
   else
       warndlg('The Excel table has different dimension that files readed')
       delete(handles.hwait)
   end

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%Calculate
%%%%%%%%%%%%%%%%%
% --- Executes on button press in calculatepush.
function calculatepush_Callback(hObject, eventdata, handles)

[Control]=ControlVar(handles);

if Control==size(get(handles.sedimentaTable, 'data'),1)*3;
    %Read input data
    if handles.excelimport==1;%desde excel
        tableData = get(handles.sedimentaTable, 'data');
        Df=cell2mat(tableData(:,1))./1000;
        Csf=cell2mat(tableData(:,2));
        Dc=cell2mat(tableData(:,3))./1000;
    elseif handles.excelimport==0;%input manual
        tableData = get(handles.sedimentaTable, 'data');
        Df=str2double(tableData(:,1))./1000;
        Csf=str2double(tableData(:,2));
        Dc=str2double(tableData(:,3))./1000;
    end

    %Extrapolation data
    CssExtraM=interpoconcentra(handles);
    VelExtraM=interpoveloci(handles);
    LatExtra=interpolateral(handles);
    
    %ADCP data
    ADCPdata=getappdata(0, 'ADCPdata');
    guiparams = getappdata(0,'guiparams');
    
for i=1:handles.numfile
    
    % Get the Application data:
    % -------------------------

    Inst{i}= guiparams{i}.Inst;
    Sensor{i}= guiparams{i}.Sensor;
    Cfg{i}= guiparams{i}.Cfg;
    V   = guiparams{i}.V;
    A   = guiparams{i}.A;
    X   = guiparams{i}.xutm;
    Y   = guiparams{i}.yutm;
    Filename = guiparams{i}.FileName;


    waitmessage=['Processing data...' Filename];
    handles.hwait = waitbar(0,waitmessage,'Name','ASET');
    guidata(hObject,handles) 
    setappdata(handles.hwait,'canceling',0)

    S.Df=Df(i);%coarse diameter
    S.Csf=Csf(i);%fine concetration
    S.Dc=Dc(i);%coarse diameter
    S.Css=10;%coarse concetration default value

    M=1;%Calculate modules (single and multiselect)
    [St{i}]=ASET_ParamStatic(M,V);
    Mont=0;%without MonteCarlo simulation
    iterations=[];%MonteCarlo iterations
    R=[];%data of Calibration
    [P{i}]=ASET_Parameters(M,V,S,R,Inst{i},Cfg{i},Sensor{i},St{i},ADCPdata,Mont,1);

    %%
    waitbar(0.1,handles.hwait);
    %%
    [Cut_vel{i}]=ASET_ArrayResizeVel(V,St{i});
    [Cut_back{i}]=ASET_ArrayResizeBack(V,St{i});
    handles.C_vel{i}=Cut_vel{i};
    handles.C_back{i}=Cut_back{i};
    guidata(hObject,handles)
    waitbar(0.4,handles.hwait);

    %%

    [Css{i},Scb{i}]=ASET_Concentration(M,V,S,Cut_back{i},R,St{i},P{i});
    waitbar(0.5,handles.hwait);

    %%

    waitbar(0.6,handles.hwait);
%%
    if (strcmp(CssExtraM{i},'NONE')) & (strcmp(VelExtraM{i},'NONE')) | (strcmp(VelExtraM{i},'NONE')) 
        VelExtra{i} = [];
        C_vel{i} = []; 
        if (strcmp(VelExtraM{i},'NONE')) & (strcmp(CssExtraM{i},'NONE'))
            C_back{i} = [];
            CssExtra{i} = [];
        elseif (strcmp(VelExtraM{i},'NONE'))       
            [C_back{i}]=ASET_ArrayResizeBackExtrap(V,St{i},Cut_back{i});
            [CssExtra{i}]=ASET_ExtrapCss(V,S,Css{i},Cut_back{i},St{i},CssExtraM{i},P{i},C_back{i});
        end
        [Dis{i}]=ASET_DischargeLiquid(V,VelExtra{i},C_vel{i},St{i},A,guiparams{i}.readformatfile);
        [Meas{i}]=ASET_TransportMeas(Css{i},V,S,St{i},Dis{i});
        [NoMeas{i}]=ASET_TransportNoMeas(V,S,CssExtra{i},VelExtra{i},C_back{i},C_vel{i},St{i},Dis{i},Cut_vel{i},Cut_back{i});
%        NoMeas{i}=[];  
    elseif (strcmp(CssExtraM{i},'NONE'))%without concetration extrapolation
        [C_vel{i}]=ASET_ArrayResizeVelExtrap(V,St{i},Cut_vel{i});
        [VelExtra{i}]=ASET_ExtrapV(V,Cut_vel{i},St{i},VelExtraM{i},C_vel{i});
        [C_back{i}]=[];
        [CssExtra{i}]=[];
        [Dis{i}]=ASET_DischargeLiquid(V,VelExtra{i},C_vel{i},St{i},A,guiparams{i}.readformatfile);
        [Meas{i}]=ASET_TransportMeas(Css{i},V,S,St{i},Dis{i});
        [NoMeas{i}]=ASET_TransportNoMeas(V,S,CssExtra{i},VelExtra{i},C_back{i},C_vel{i},St{i},Dis{i},Cut_vel{i},Cut_back{i});
        %NoMeas{i}=[];    
    else %with extrapolations both variables
        [C_vel{i}]=ASET_ArrayResizeVelExtrap(V,St{i},Cut_vel{i});
        [VelExtra{i}]=ASET_ExtrapV(V,Cut_vel{i},St{i},VelExtraM{i},C_vel{i});
        [C_back{i}]=ASET_ArrayResizeBackExtrap(V,St{i},Cut_back{i});
        [CssExtra{i}]=ASET_ExtrapCss(V,S,Css{i},Cut_back{i},St{i},CssExtraM{i},P{i},C_back{i});
        [Dis{i}]=ASET_DischargeLiquid(V,VelExtra{i},C_vel{i},St{i},A,guiparams{i}.readformatfile);
        [Meas{i}]=ASET_TransportMeas(Css{i},V,S,St{i},Dis{i});
        [NoMeas{i}]=ASET_TransportNoMeas(V,S,CssExtra{i},VelExtra{i},C_back{i},C_vel{i},St{i},Dis{i},Cut_vel{i},Cut_back{i});
    end
    
    %
    waitbar(0.6,handles.hwait);
    
    %Concanate all data
    [Array{i}]=ASET_Concatenate(V,Css{i},Cut_vel{i},Cut_back{i},VelExtra{i},CssExtra{i},St{i},Meas{i},NoMeas{i},C_vel{i},C_back{i});
     
    %Lateral Extrapolation. Only possible when the input file is ASCII
    if LatExtra{i}==1
        QEdge{i}.Gss.left=0;
        QEdge{i}.Gss.right=0;
        QEdge{i}.Gw.left=0;
        QEdge{i}.Gw.right=0;
        QEdge{i}.Q.right=0;
        QEdge{i}.Q.left=0;
    else %Extrapolation and ASCI file loaded
        if LatExtra{i}==2
            Method.Edge=0;%triangular method
        elseif LatExtra{i}==3
            Method.Edge=1;%rectangular method
        end
        %Find file and determinate left and right distance
        Edge.leftdistance = guiparams{i}.Edge.leftdistance;%
        Edge.rightdistance = guiparams{i}.Edge.rightdistance;%
        [QEdge{i}]=ASET_LateralExtraDis(Method,Array{i},V,S,Edge);
    end
       
    %Final calculate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Coarse material transport
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NoMeas{i}.Surf.GssT=nansum(nansum(NoMeas{i}.Surf.coarse));%Surface zone transport
    Meas{i}.GssT=nansum(nansum(Meas{i}.coarse));%Measurement zone trasport
    NoMeas{i}.Bottom.GssT=nansum(nansum(NoMeas{i}.Bottom.coarse));%Bottom zone trasport

    Array{i}.GssT=QEdge{i}.Gss.left+QEdge{i}.Gss.right+NoMeas{i}.Surf.GssT+...
        Meas{i}.GssT+NoMeas{i}.Bottom.GssT;%Total transport

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fine material transport
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
    NoMeas{i}.Surf.GwT=nansum(nansum(NoMeas{i}.Surf.fine));%Surface zone transport
    Meas{i}.GwT=nansum(nansum(Meas{i}.fine));%Measurement zone trasport
    NoMeas{i}.Bottom.GwT=nansum(nansum(NoMeas{i}.Bottom.fine));%Bottom zone trasport

    Array{i}.GwT=NoMeas{i}.Surf.GwT+Meas{i}.GwT+NoMeas{i}.Bottom.GwT+...
        QEdge{i}.Gw.left+QEdge{i}.Gw.right;%Washload

    %Write results    
    FinalTable(i,:)=[round(NoMeas{i}.Surf.GssT,2) round(Meas{i}.GssT,2) ...
        round(NoMeas{i}.Bottom.GssT,2) round(Array{i}.GssT,2)...
        round(Array{i}.GwT,2) round(Dis{i}.Total+QEdge{i}.Q.right+QEdge{i}.Q.left,2)];

    guidata(hObject,handles)


    %Storage data
    guiparams{i}.S=S;
    guiparams{i}.V=V;
    guiparams{i}.Dis=Dis{i};
    guiparams{i}.VelExtraM=VelExtraM(i);
    guiparams{i}.CssExtraM=CssExtraM(i);
    guiparams{i}.LatExtra=LatExtra(i);
    guiparams{i}.St=St{i};
    guiparams{i}.P=P{i};
    guiparams{i}.C_vel=C_vel{i};
    guiparams{i}.C_back=C_back{i};
    guiparams{i}.Css=Css{i};
    guiparams{i}.Meas=Meas{i};
    guiparams{i}.Array=Array{i};
    guiparams{i}.QEdge=QEdge{i};
    guiparams{i}.NoMeas=NoMeas{i};
    guiparams{i}.VelExtra=VelExtra{i};
    guiparams{i}.CssExtra=CssExtra{i};

    clear R V A X Y    

    waitbar(1,handles.hwait)

    delete(handles.hwait)

end

%Storage data
setappdata(0, 'guiparams', guiparams)

set(handles.results,'Data',FinalTable)

%Show the results 
showSummaryTable(guiparams,handles)

% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'Calculate Complete '};
statusLogging(handles.LogWindow, log_text)


clear S Var CssExtra VelExtraM CssExtraM St P C_vel C_back Css Meas NoMeas...
    VelExtra Array LatExtra 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%Extra Functions
%%%%%%%%%%%%%%%%%

function set_enable(handles,enable_state)
switch enable_state
    case 'init'
        set(handles.results,'Enable','off');
        set(handles.inputfiles,'Enable','off');
        set(handles.sedimentaTable,'Enable','off');
        set(handles.calculatepush,'Enable','off');
        set(handles.adcpdataoption,'Enable','off');
        set(handles.importdata,'Enable','off');
        set(handles.quitfiles,'Enable','off');
        set(handles.exportfunction,'Enable','off')
        set(handles.SumT,'Enable', 'off')
        set(handles.sedimentaTable, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) 6]));
        set(handles.results, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) 5]));
        set(handles.inputfiles, 'Data', cell([size(get(handles.inputfiles,'Data'),1) 1]));
        set(handles.SumT, 'Data', cell([size(get(handles.sedimentaTable,'Data'),1) 11]));
    case 'loadfiles'
        set(handles.results,'Enable','on');
        set(handles.inputfiles,'Enable','on');
        set(handles.sedimentaTable,'Enable','on');
        set(handles.calculatepush,'Enable','on');
        set(handles.adcpdataoption,'Enable','on');
        set(handles.importdata,'Enable','on');
        set(handles.exportfunction,'Enable','on')
        set(handles.SumT,'Enable','on')
        set(handles.quitfiles,'Enable','on');
    otherwise
end


function CssExtraM=interpoconcentra(handles)
%Read extrapolation concentration selected and modify
tableData = get(handles.sedimentaTable, 'data');  
handles.CssExtraM=cell(handles.numfile,1);

for i=1:handles.numfile
    if strcmp(tableData(i,5),'None')
            CssExtraM{i,1}='NONE';%Constant
    elseif strcmp(tableData(i,5),'Constant')
            CssExtraM{i,1}='CON';%Linear Method
    elseif strcmp(tableData(i,5),'3Pt. Slope')
            CssExtraM{i,1}='LSI';%Linear Method
    elseif strcmp(tableData(i,5),'Rouse')
            CssExtraM{i,1}='RSI';%Rouse Method
    end

end

function VelExtraM=interpoveloci(handles)
%Read extrapolation velocity selected and modify
tableData = get(handles.sedimentaTable, 'data');  
handles.VelExtraM=cell(handles.numfile,1);

for i=1:handles.numfile
    if strcmp(tableData(i,4),'None')
            VelExtraM{i,1}='NONE';
    elseif strcmp(tableData(i,4),'Constant')
            VelExtraM{i,1}='CON';%Constant 
    elseif strcmp(tableData(i,4),'3Pt. Slope')
            VelExtraM{i,1}='LVI';%3Pt Slope 
    elseif strcmp(tableData(i,4),'Law of the Wall')
            VelExtraM{i,1}='LPVI';%Law of the wall method
    end
end

function LatExtra=interpolateral(handles)
 %Extrapolation methods lateral
for i=1:handles.numfile
    tableData = get(handles.sedimentaTable, 'data');  
    handles.LatExtra=cell(handles.numfile,1);
    if strcmp(tableData(i,6),'None')
            LatExtra{i,1}=1;
    elseif strcmp(tableData(i,6),'Triangular')
            LatExtra{i,1}=2;% Triangular Method 
    elseif strcmp(tableData(i,6),'Rectangular')
            LatExtra{i,1}=3;% Rectangular Method
    end
end


% --------------------------------------------------------------------
function closefunction_Callback(hObject, eventdata, handles)
close


% --- Executes on selection change in LogWindow.
function LogWindow_Callback(hObject, eventdata, handles)
% Empty


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

function [Control]=ControlVar(handles)

for i=1:size(get(handles.sedimentaTable, 'data'),1);
    tableData=get(handles.sedimentaTable, 'data');
    if isempty(cell2mat(tableData(i,1)))
        wrong_data=warndlg('Please complete the Fine Sediment Diameter box',...
        'WARNING'); 
        Control1(i)=0;
    else
        Control1(i)=1;
    end
    if isempty(cell2mat(tableData(i,2)))
        wrong_data=warndlg('Please complete the Coarse Sediment Diameter box',...
        'WARNING');
        Control2(i)=0;
    else
        Control2(i)=1;
    end
    if isempty(cell2mat(tableData(i,3)))
        wrong_data=warndlg('Please complete the Fine Sediment Concentration box',...
        'WARNING'); 
        Control3(i)=0;
    else
        Control3(i)=1;
    end
    
end

    Control=nansum(Control1)+nansum(Control2)+nansum(Control3);
    
    
function showSummaryTable(guiparams,handles)
% This function create the array for the summary table in mulstiselect GUI
for t=1:handles.numfile
SummaryT(t,:)=[guiparams{t}.V.navref sprintf('%.2f',guiparams{t}.V.mcsDist(1,end))...
    sprintf('%.2f',nanmean(guiparams{t}.Array.Bed)) sprintf('%.2f',nanmean(nanmean(guiparams{t}.Array.uMag))) ...
    sprintf('%.2f',guiparams{t}.Dis.surfT) sprintf('%.2f',guiparams{t}.Dis.bottomT)...
    sprintf('%.2f',guiparams{t}.QEdge.Q.left) sprintf('%.2f',guiparams{t}.QEdge.Q.right) ...
    sprintf('%.2f',nanmean(nanmean(guiparams{t}.Array.Css))*1000)...
    sprintf('%.2f',guiparams{t}.QEdge.Gss.left) sprintf('%.2f',guiparams{t}.QEdge.Gss.right)];
end
set(handles.SumT,'Data',SummaryT)

% --------------------------------------------------------------------
function set_Callback(hObject, eventdata, handles)
% empty


% --------------------------------------------------------------------
function unitsset_Callback(hObject, eventdata, handles)
% empty


% --------------------------------------------------------------------
function unitsMetric_Callback(hObject, eventdata, handles)
handles.units=0;
guidata(hObject,handles)

% Update the GUI:
% ---------------
set(handles.metricunits, 'Checked','on')
set(handles.englishunits,'Checked','off')

% --------------------------------------------------------------------
function unitsEnglish_Callback(hObject, eventdata, handles)
handles.units=1;
guidata(hObject,handles)

% Update the GUI:
% ---------------
set(handles.metricunits, 'Checked','off')
set(handles.englishunits,'Checked','on')


% --- Executes when selected cell(s) is changed in inputfiles.
function inputfiles_CellSelectionCallback(hObject, eventdata, handles)
%define the number of file selected
handles.quitfile=eventdata.Indices(:,1);
guidata(hObject,handles)
