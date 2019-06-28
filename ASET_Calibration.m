function varargout = ASET_Calibration(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function Calibrate the signal using the methodlogt propused by 
% Szupiany et al (2018). The parameter determinated by this function is the
% kcEr. 

%by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dominguez Ruben 5/10/2018

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_Calibration_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_Calibration_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%
%Open Function
%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before ASET_Calibration is made visible.
function ASET_Calibration_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ASET_Calibration
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

%Control data
handles.controlread=0;%Metodo de lectura
count=0;
setappdata(0,'count',count)
handles.readadcpdata=0;%read adcp data
handles.SAMPLE=0;%Single select
handles.excelimport=0;
handles.LoadPD0=0;
handles.first=1;%firs read
handles.quitfile=0;
guidata(hObject, handles);
set(handles.inputfilename, 'Data',cell(4,1));


set_enable(handles,'init')%set initial configuration

%Initial Panel
set(handles.figure1,'Name',['ASET: Calibration'], ...
    'DockControls','off')
% UIWAIT makes ASET_Calibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ASET_Calibration_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%
%Toolbar Menu
%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% empty
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function newproject_Callback(hObject, eventdata, handles)
%Delete all data for a new project

handles.first=1;
guidata(hObject,handles)
set_enable(handles,'init')
set(handles.addfiles,'Enable','on');
% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%----------- ' datestr(now) ' ------------%'];...
    'New Project'};
statusLogging(handles.LogWindow, log_text)

function closefunction_Callback(hObject, eventdata, handles)
close

%%%%%%%%%%%%%%%%%%%%%%%%
%Main Panel
%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in addfiles.
function addfiles_Callback(hObject, eventdata, handles)


if handles.first==1;%first loaded
     %Read .PD0 file
    [filenamepd0,pathnamepd0] = uigetfile({'*.PD0';'*.000';'*.*'},'Select .PD0 files','MultiSelect','on');
    
    if iscell(filenamepd0)
        if iscell(filenamepd0)
            handles.length=size(filenamepd0,2);
        else
            handles.length=size(filenamepd0,1);
        end
        
        if handles.length==1
            infile=fullfile(pathnamepd0,filenamepd0);
            set(handles.inputfilename,'Data',filenamepd0');
            handles.infile=infile;
            handles.filenamepd0=filenamepd0;
            guidata(hObject,handles)
        else
            for y=1:handles.length
                infile(y)=fullfile(pathnamepd0,filenamepd0(y));
            end
            set(handles.inputfilename,'Data',filenamepd0');
            handles.infile=infile;
            handles.filenamepd0=filenamepd0;
            guidata(hObject,handles)
        end
        
        handles.Calibration=1;%Calibration to adcpdata script
        handles.draft=0;%Without draft data
        handles.R.multi=1;%multiselect on
        handles.first=0;
        handles.pathPD0=pathnamepd0;
        
        handles.noenablekcEr=1;%deactive kcEr box
        guidata(hObject,handles)

        ASET_ADCPData(handles)

        handles.first=0;
        
        set_enable(handles,'loadfiles')
          
        clear infile
    end
elseif handles.first==0;%second loaded

     if iscell(handles.infile)
        infile=handles.infile;
        filenamepd0pre=handles.filenamepd0;
     else
        infile={handles.infile};
        filenamepd0pre={handles.filenamepd0};
     end
    
    %Read .PD0 file
    [filenamepd0,pathnamepd0] = uigetfile({'*.PD0';'*.000';'*.*'},'Select .PD0 files',handles.pathPD0,'MultiSelect','on');
    
    if iscell(filenamepd0)
        handles.length=size(filenamepd0,2);
    else
        handles.length=size(filenamepd0,1);
    end
    
        if handles.length==1
            infile1=fullfile(pathnamepd0,filenamepd0);
            
            if length(filenamepd0pre)==1
                infile=[infile';infile1];
                handles.filenamepd0=[handles.filenamepd0;{filenamepd0}];
            else
                infile=[infile';infile1];
                handles.filenamepd0=[handles.filenamepd0';{filenamepd0}];
            end
            set(handles.inputfilename,'Data',handles.filenamepd0)        
            guidata(hObject,handles)
        else
        
        for y=1:handles.length
            infile1(y)=fullfile(pathnamepd0,filenamepd0(y)); 
        end
            
            if length(filenamepd0pre)==1
                handles.filenamepd0=[str2mat(handles.filenamepd0);filenamepd0'];
                infile=[infile';infile1'];
            else
                handles.filenamepd0=[handles.filenamepd0';filenamepd0'];
                infile=[infile';infile1'];
            end
            set(handles.inputfilename,'Data',handles.filenamepd0);

            guidata(hObject,handles)
        end      
        
        handles.Calibration=1;%Calibration on to adcpdata script
        handles.draft=0;%Without draft data
        handles.R.multi=1;%multiselect on

        handles.noenablekcEr=1;%deactive kcEr box
        handles.first=0;
        handles.Calibration=1;
        handles.length=length(handles.filenamepd0);
        handles.infile=infile;
        guidata(hObject,handles)
        
        ASET_ADCPData(handles)
        set_enable(handles,'loadfiles')  
        clear infile

end


% --------------------------------------------------------------------
function readpd0_Callback(hObject, eventdata, handles)
%empty
 
%%%%%%%%%%%%%%%%%% 
%Load PD0 files
%%%%%%%%%%%%%%%%%%
function loadfilespd0_Callback(hObject, eventdata, handles)
%Load data   
set_enable(handles,'loadpd0')
infile=handles.infile;

filenamepd0=handles.filenamepd0;

M=2;%Calibration Module

ADCPdata=getappdata(0, 'ADCPdata');
%  
 if handles.length>1%more than one file
     for i=1:handles.length
        [R{i},Inst{i},Cfg{i},Sensor{i}]=ASET_ReadPd0(infile{i},ADCPdata.R.draft(i),M); 

        % Push messages to Log Window:
        % ----------------------------
        log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'PD0 File Loaded';[cell2mat(infile(i))]};
        statusLogging(handles.LogWindow, log_text)
     end
 else %one file
    [R,Inst,Cfg,Sensor]=ASET_ReadPd0(infile,ADCPdata.R.draft,M); 

    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'PD0 File Loaded';[cell2mat(infile)]};
    statusLogging(handles.LogWindow, log_text)
 end
            
%Store pd0 data
setappdata(0,'R',R)
setappdata(0,'Inst',Inst)
setappdata(0,'Cfg',Cfg)
setappdata(0,'Sensor',Sensor)

%Built the table sedimentology
celltable=cell(handles.length,7);
celltable(:,3:7)={''};
celltable(:,2)={'1'};
celltable(:,1)={filenamepd0{:}};
handles.SAMPLE=str2double(celltable(:,2));

set(handles.uitablesedimentology,'Data',celltable)

if handles.length>1%more than two files
     handles.Inst=Inst{1};
     handles.R.draft=R{1}.draft;
else% one file
     handles.Inst{1}=Inst;
     handles.R.draft=R.draft;
end
  
if size(celltable,1)<3
    wrong_data=warndlg('You need add at least 3 samples to calibration',...
    'WARNING'); 
    set(handles.uitablesedimentology, 'Data',...
        cell(size(get(handles.uitablesedimentology,'Data'))));
else

end
  
%Change the draft and rssi data
m=1;
for i=1:handles.length
    for t=1:handles.SAMPLE(i)
        ADCPdata1.rssi_beam1(m)=ADCPdata.rssi_beam1(i);
        ADCPdata1.rssi_beam2(m)=ADCPdata.rssi_beam2(i);
        ADCPdata1.rssi_beam3(m)=ADCPdata.rssi_beam3(i);
        ADCPdata1.rssi_beam4(m)=ADCPdata.rssi_beam4(i);
        ADCPdata1.R.draft(m)=ADCPdata.R.draft(i);
        m=m+1;
    end
end

%Re-writer variables
ADCPdata.rssi_beam1=ADCPdata1.rssi_beam1;
ADCPdata.rssi_beam2=ADCPdata1.rssi_beam2;
ADCPdata.rssi_beam3=ADCPdata1.rssi_beam3;
ADCPdata.rssi_beam4=ADCPdata1.rssi_beam4;

set(handles.addfiles,'Enable','off');
guidata(hObject,handles)      
setappdata(0,'ADCPdata',ADCPdata)  
%Control complete calculate
handles.LoadPD0=1;
guidata(hObject, handles);
   
 
% --------------------------------------------------------------------
function savefunction_Callback(hObject, eventdata, handles)
%Empty


 % --------------------------------------------------------------------
function matfileexport_Callback(hObject, eventdata, handles)
%Save the calibration file
Cal=getappdata(0,'Cal');
S=getappdata(0,'S');
P=getappdata(0,'P');
Graph=getappdata(0,'Graph');
R=getappdata(0,'R');
        
[file,path] = uiputfile('*.mat','Save file');
if file==0
else
    save([path file], 'Cal','S','P','Graph','R');

% Push messages to Log Window:
    % ----------------------------
    log_text = {...
        '';...
        ['%--- ' datestr(now) ' ---%'];...
        'Complete Export .mat file' ;[path file]};
    statusLogging(handles.LogWindow, log_text)
end

% --------------------------------------------------------------------
function excelfilesexport_Callback(hObject, eventdata, handles)
%%%%
Cal=getappdata(0,'Cal');
S=getappdata(0,'S');
Graph=getappdata(0,'Graph');
P=getappdata(0,'P');
[file,path] = uiputfile('*.xlsx','Save *.xlsx file');

if file==0
else
    outfile = fullfile(path,file);  
    hwait = waitbar(0,'Exporting Excel File...','Name','ASET');

    if P{1}.Ac.FREQ==1228.8
        FREQ=1200;
    elseif P{1}.Ac.FREQ==614
        FREQ=1200;
    end

    sout = {...
            'ASET: Summary Calibration Data' '' '' '' '' '';...
            'Date Processed: ' datestr(now) '' '' '' '';...
            'Path and File' '' '' '' '' outfile;
            '' '' '' '' '' '';...
            'Calibration Analisys:' '' '' '' '' '';...
            'Numbers of Samplers:' '' '' '' '' length(P);...
            'ADCP Frenquency [kHz]:' '' '' '' '' FREQ;...
            'Draft [m]:' '' '' '' '' P{1}.Ac.Draft_ADCP;...
            'kcEr:' '' '' '' '' Cal.kcEr_mean;...
            'RSSI Slope Beam 1' '' '' '' '' P{1}.Ac.beam_1;...
            'RSSI Slope Beam 2' '' '' '' '' P{1}.Ac.beam_2;...
            'RSSI Slope Beam 3' '' '' '' '' P{1}.Ac.beam_3;...
            'RSSI Slope Beam 4' '' '' '' '' P{1}.Ac.beam_4;...
            '' '' '' '' '' '';...
            '' '' '' '' '' '';...
            'Fitting Data' '' '' '' '' '';...
            'ST vs Ms2 Fitting' '' '' Graph.coef1(1) 'x' Graph.coef1(2);...
            'Ms2_M vs Ms2_C Fitting' '' '' Graph.coef2(1) 'x' Graph.coef2(2);...
            'R-Squared ST vs Ms2' '' '' '' '' Graph.r1;...
            'R-Squared Ms2_M vs Ms2_C' '' '' '' '' Graph.r2;...
            };
        xlswrite(outfile,sout,'ASET Calibration Summary','A1');
        waitbar(1/3,hwait)

    %Find the data
    for i=1:length(P)
      if isnan(Cal.ST{i}(Cal.indice(i))) | isnan(Cal.kcEr{i}(Cal.indice(i))) | isnan(Cal.SCB{i}(Cal.indice(i))) | isnan(Cal.Css{i}(Cal.indice(i))) 
        ST(i,1)=-999;
        kcEr(i,1)=-999;
        SCB(i,1)=-999;
        Css(i,1)=-999;
      else
        ST(i,1)=Cal.ST{i}(Cal.indice(i));
        kcEr(i,1)=Cal.kcEr{i}(Cal.indice(i));
        SCB(i,1)=Cal.SCB{i}(Cal.indice(i));
        Css(i,1)=Cal.Css{i}(Cal.indice(i));
      end
    end

    TABLEData=get(handles.uitablesedimentology,'Data');

    s1=[TABLEData num2cell([ST Css kcEr SCB])];

    s1headers = {'Filename' 'Nº Sample' 'Fine Diameter [mm]' 'Ms1 Concetration [mg/l]'...
    'Coarse Diameter [mm]' 'Ms2 Concentration Measured  [mg/l]' 'Depth from surface [m]' ...
    'ST [dB]' 'Ms2 Calculate [mg/l]' 'kcEr' 'SCB'};

    pvout1 = vertcat(s1headers,s1);
    xlswrite(outfile,pvout1, 'ASET Calibration Data');
    waitbar(3/3,hwait)
    delete(hwait)

    % Push messages to Log Window:
    % ----------------------------
    log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'Complete Export .xlsx file';outfile};
    statusLogging(handles.LogWindow, log_text)
end


% --------------------------------------------------------------------
function figuresfunction_Callback(hObject, eventdata, handles)
%Empty


% --------------------------------------------------------------------
function exportgraph1_Callback(hObject, eventdata, handles)
Graph=getappdata(0,'Graph');

h=figure('Name','ASET: ST vs Ms2','NumberTitle','off');
set(h, 'Position', [10 10 800 400])
AxesH = handles.calgraph(1);
fig=copyobj(AxesH,h);
set(fig,'Units', 'normalized', 'Position', [.1 .11 .85 .8]);

%show equation    
eq1 = uicontrol('style','text',...
'position',[100 340 100 20],...
'fontsize',12);
set(eq1,'string',['y=',num2str(round(Graph.coef1(1)*100)/100),'x',num2str(round(Graph.coef1(2)*100)/100)]);


% --------------------------------------------------------------------
function exportgraph2_Callback(hObject, eventdata, handles)
Graph=getappdata(0,'Graph');

h=figure('Name','ASET: Ms2_M vs Ms2_C','NumberTitle','off');
set(h, 'Position', [10 10 800 400])
AxesH = handles.calgraph2(1);
fig=copyobj(AxesH,h);
set(fig,'Units', 'normalized', 'Position', [.1 .11 .85 .8]);

%show equation    
eq1 = uicontrol('style','text',...
'position',[100 340 100 20],...
'fontsize',12);
set(eq1,'string',['y=',num2str(round(Graph.coef2(1)*100)/100),'x',num2str(round(Graph.coef2(2)*100)/100)]);


%%%%%%%%%%%%%%%%%%%%%%%%%
%Sedimentology Table Data
%%%%%%%%%%%%%%%%%%%%%%%%%

function finediameter_Callback(hObject, eventdata, handles)
%Empty


% --- Executes during object creation, after setting all properties.
function finediameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to finediameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fineconcentration_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function fineconcentration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fineconcentration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function coursediameter_Callback(hObject, eventdata, handles)
% Empty

% --- Executes during object creation, after setting all properties.
function coursediameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coursediameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function courseconcentration_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function courseconcentration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to courseconcentration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD0 Load
%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on key press with focus on loadfilespd0 and none of its controls.
function loadfilespd0_KeyPressFcn(hObject, eventdata, handles)
% Empty


% --- Executes on button press in importsed.
function importsed_Callback(hObject, eventdata, handles)
%Import Sedimentology data from Excel files
[filenamexls,pathnamexls] = uigetfile({'*.xlsx';'*.*'},'Select *.xlsx file');

if filenamexls==0
else
    hwait = waitbar(0,'Importing Excel File...','Name','ASET');
     waitbar(1/3,hwait)

    if ischar(filenamexls)
       infile=[pathnamexls filenamexls];
       xls=num2cell(xlsread(infile));
       if ischar(xls{1,1})
           xls{1,:}=[];
       else
       end
    end
    
    if size(xls,1)==handles.length
        
        waitbar(2/3,hwait)

        datapd0=get(handles.uitablesedimentology, 'data');
        filename=datapd0(:,1);

         set(handles.uitablesedimentology, 'data',[filename num2cell(handles.SAMPLE) xls]);
         handles.excelimport=1;
         guidata(hObject, handles);

        waitbar(3/3,hwait)
        delete(hwait)
    else
        warndlg('The Excel table has different dimension that files readed')
        delete(hwait)
    end
end




%%%%%%%%%%%%%%%%%%%%%%%
% Graph data
%%%%%%%%%%%%%%%%%%%%%%%

function r1_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function r1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r2_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function r2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Press Calibration Buttton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in calibrate.
function calibrate_Callback(hObject, eventdata, handles)

 hwait = waitbar(0,'Calibration ...','Name','ASET');
 
%Read input data
if handles.excelimport==1;%desde excel
    tableData = get(handles.uitablesedimentology, 'data');
    S.Df=cell2mat(tableData(:,3))./1000;
    S.Csf=cell2mat(tableData(:,4));
    S.Dc=cell2mat(tableData(:,5))./1000;
    S.Css=cell2mat(tableData(:,6));
    S.depthdata=cell2mat(tableData(:,7));
elseif handles.excelimport==0;%input manual
    tableData = get(handles.uitablesedimentology, 'data');
    S.Df=str2double(tableData(:,3))./1000;
    S.Csf=str2double(tableData(:,4));
    S.Dc=str2double(tableData(:,5))./1000;
    S.Css=str2double(tableData(:,6));
    S.depthdata=str2double(tableData(:,7));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Method 2
M=2; %Calibration module
V=[];%No Data
C=[]; %No Data

R1 = getappdata(0,'R'); 
Inst1 = getappdata(0,'Inst');
Cfg1 = getappdata(0,'Cfg');
Sensor1 = getappdata(0,'Sensor');
%ADCP data
ADCPdata=getappdata(0, 'ADCPdata');
p=1;
    
for u=1:length(R1)
     if handles.SAMPLE(u)>1 %more than one sample per each PD0 File  
         for i=1:length(handles.SAMPLE(u))
                R{p}.draft=R1{i}.draft;
                R{p}.Depth_m=R1{i}.Depth_m;
                R{p}.Back.Beam1(:,:,1)=R1{i}.Back.Beam1(:,:,1);
                R{p}.Back.Beam2(:,:,1)=R1{i}.Back.Beam2(:,:,1);
                R{p}.Back.Beam3(:,:,1)=R1{i}.Back.Beam3(:,:,1);
                R{p}.Back.Beam4(:,:,1)=R1{i}.Back.Beam4(:,:,1);
                Inst(p)= Inst1(i);
                Cfg(p)= Cfg1(i);
                Sensor(p)= Sensor1(i);
                p=p+1;
         end
     elseif handles.SAMPLE(u)==1;%one sample one PD0 File
            R{p}.draft=R1{u}.draft;
            R{p}.Depth_m=R1{u}.Depth_m;
            R{p}.Back.Beam1(:,:,1)=R1{u}.Back.Beam1(:,:,1);
            R{p}.Back.Beam2(:,:,1)=R1{u}.Back.Beam2(:,:,1);
            R{p}.Back.Beam3(:,:,1)=R1{u}.Back.Beam3(:,:,1);
            R{p}.Back.Beam4(:,:,1)=R1{u}.Back.Beam4(:,:,1);
            Inst{p}= Inst1{u};
            Cfg{p}= Cfg1{u};
            Sensor{p}= Sensor1{u};
            p=p+1;
     end
end
 

clear R1 Inst1 Cfg1 Sensor1
    
for t=1:size(R,2)
    [St{t}]=ASET_ParamStatic(M,V,R{t});
    Mont=0;%Without Monte Carlo
    iterations=[];%MonteCarlo iterations

    [P{t}]=ASET_Parameters(M,V,S,R,Inst{t},Cfg{t},Sensor{t},St{t},ADCPdata,Mont,t);
end

waitbar(0.25,hwait)

%Init Calibration
[~,~,Cal]=ASET_Concentration(M,V,S,C,R,St,P);
%End calibration

%Graph results
[Graph]=ASET_GraphCal(S,Cal,R,P);


waitbar(0.5,hwait)

%Store Data
setappdata(0,'Cal',Cal)
setappdata(0,'S',S)
setappdata(0,'P',P)
setappdata(0,'Graph',Graph)

%Graphics builder
%clear graph
cla(handles.calgraph)
cla(handles.calgraph2)
waitbar(0.75,hwait)
%Picture 1
axes(handles.calgraph)
plot(Graph.x1,Graph.y1,'or')
hold on
plot(Graph.x1pred,Graph.y1pred,'-b','LineWidth',1)
xlabel('ST [dB]');ylabel('log(Ms2_M*ks)')
hold off

if Graph.r1<0.7
    set(handles.r1,'String',num2str(round(Graph.r1*100)/100));
    set(handles.r1,'BackgroundColor','red');
else
    set(handles.r1,'String',num2str(round(Graph.r1*100)/100));
    set(handles.r1,'BackgroundColor','white');
end
    
%show equation    
posax1=get(handles.calgraph,'position');
eq1 = uicontrol('style','text',...
'position',[posax1(1)+300 posax1(2)+600 posax1(3)+50 posax1(4)],...
'fontsize',12);
if Graph.coef1(2)>0
    set(eq1,'string',['y=',num2str(round(Graph.coef1(1)*100)/100),'x+',num2str(round(Graph.coef1(2)*100)/100)]);
elseif Graph.coef1(2)<0
    set(eq1,'string',['y=',num2str(round(Graph.coef1(1)*100)/100),'x',num2str(round(Graph.coef1(2)*100)/100)]);
end

%Picture 2
axes(handles.calgraph2)
plot(Graph.x2,Graph.y2,'ob')
hold on
plot(Graph.x2pred,Graph.y2pred,'-r','LineWidth',1)
xlabel('10*log(Ms2_M)');ylabel('10*log(Ms2_C)')
hold off

if Graph.r2<0.7
    set(handles.r2,'String',num2str(round(Graph.r2*100)/100));
    set(handles.r2,'BackgroundColor','red');
else
    set(handles.r2,'String',num2str(round(Graph.r2*100)/100));
    set(handles.r2,'BackgroundColor','white');
end

%show equation
posax2=get(handles.calgraph2,'position');
eq2 = uicontrol('style','text',...
'position',[posax2(1)+600 posax2(2)+600 posax2(3)+50 posax2(4)],...
'fontsize',12);

if Graph.coef2(2)>0
    set(eq2,'string',['y=',num2str(round(Graph.coef2(1)*100)/100),'x+',num2str(round(Graph.coef2(2)*100)/100)]);
elseif Graph.coef2(2)<0
    set(eq2,'string',['y=',num2str(round(Graph.coef2(1)*100)/100),'x',num2str(round(Graph.coef2(2)*100)/100)]);
end

set(handles.savefunction,'Enable','on');
    
%control the number of sample excluide
t=1;
h=[];
for i=1:size(R,2)
    if isnan(Cal.kcEr_puntual(i))
        h(t)=i;
        t=t+1;
    else
    end
end

waitbar(1,hwait)
delete(hwait)

% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'Complete Calculate'};
statusLogging(handles.LogWindow, log_text)
    
% Push messages to Log Window:
% ----------------------------
log_text = {...
    '';...
    ['%--- ' datestr(now) ' ---%'];...
    'kcEr';[cell2mat({Cal.kcEr_mean})] };
statusLogging(handles.LogWindow, log_text)

%Show the smples not included in the calibration (Id and name of PD0 file)
if isempty(h)
else
    for o=1:length(h)
        C{1,o}=h(o);
        C{2,o}=cell2mat(handles.filenamepd0(h(o)));
    end
    str=sprintf('\n %d: %s',C{:});
    prompt={str};
    warndlg(prompt,'Sample not included')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extra Function
%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_enable(handles,enable_state)
%Set initial and load files
switch enable_state
    case 'init'
        cla(handles.calgraph)
        cla(handles.calgraph2)
        axes(handles.calgraph2)
        grid on
        axes(handles.calgraph)
        grid on
        set(handles.uitablesedimentology,'Enable','off');
        set(handles.adcpdata,'Enable','off');
        set(handles.quitfiles,'Enable','off');
        set(handles.loadfilespd0,'Enable','off');
        set(handles.calibrate,'Enable','off');
        set(handles.importsed,'Enable','off');
        set(handles.savefunction,'Enable','off');
        set(handles.uitablesedimentology, 'Data', cell(size(get(handles.uitablesedimentology,'Data'))));
        set(handles.inputfilename, 'Data', cell(4,1));
        set(handles.inputfilename, 'columnname','Filename');
        set(handles.r2,'String','');
        set(handles.r1,'String','');
        handles.first=1;
        handles.length=[];
        handles.infile=[];
        handles.filenamepd0=[];

        posax1=get(handles.calgraph,'position');
        eq1 = uicontrol('style','text',...
        'position',[posax1(1)+300 posax1(2)+600 posax1(3)+50 posax1(4)],...
        'fontsize',12);
         set(eq1,'string','');

        posax2=get(handles.calgraph2,'position');
        eq2 = uicontrol('style','text',...
        'position',[posax2(1)+600 posax2(2)+600 posax2(3)+50 posax2(4)],...
        'fontsize',12);
        set(eq2,'string','');

        set(handles.r1,'BackgroundColor','white');
        set(handles.r2,'BackgroundColor','white');
    
    case 'loadfiles'

    set(handles.loadfilespd0,'Enable','on');
        set(handles.quitfiles,'Enable','on');
        set(handles.adcpdata,'Enable','on');

        set(handles.uitablesedimentology, 'Data', cell(size(get(handles.uitablesedimentology,'Data'))));

    case 'loadpd0'
        cla(handles.calgraph)
        cla(handles.calgraph2)
        set(handles.savefunction,'Enable','off');
        set(handles.r2,'String','');
        set(handles.r1,'String','');
        set(handles.calibrate,'Enable','on');
        set(handles.importsed,'Enable','on');
        set(handles.uitablesedimentology,'Enable','on');
        set(handles.uitablesedimentology, 'Data', cell(size(get(handles.uitablesedimentology,'Data'))));
        set(handles.r1,'BackgroundColor','white');
        set(handles.r2,'BackgroundColor','white');
    
    otherwise
end



% --- Executes on selection change in LogWindow.
function LogWindow_Callback(hObject, eventdata, handles)
% empty


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


% --- Executes on button press in adcpdata.
function adcpdata_Callback(hObject, eventdata, handles)
ADCPdata=getappdata(0, 'ADCPdata');
ASET_ADCPData(ADCPdata)


% --- Executes on button press in quitfiles.
function quitfiles_Callback(hObject, eventdata, handles)
%This funtion quit files PD0 readed

if handles.quitfile==0
   warndlg('Please select the file that you can delete')
else
    TableDATA=get(handles.inputfilename,'Data');

    TableDATA(handles.quitfile,:)=[];
    set(handles.inputfilename,'Data',TableDATA)

    if handles.LoadPD0==0;
    elseif handles.LoadPD0==1;
        TableDATASED = get(handles.uitablesedimentology,'Data');
        TableDATASED(handles.quitfile,:)=[];
        set(handles.uitablesedimentology,'Data',TableDATASED)

        R = getappdata(0,'R'); 
        Inst = getappdata(0,'Inst');
        Cfg = getappdata(0,'Cfg');
        Sensor = getappdata(0,'Sensor');

        R(handles.quitfile)=[];
        Inst(handles.quitfile)= [];
        Cfg(handles.quitfile)= [];
        Sensor(handles.quitfile)= [];

        setappdata(0,'R',R)
        setappdata(0,'Inst',Inst)
        setappdata(0,'Cfg',Cfg)
        setappdata(0,'Sensor',Sensor)
            
        handles.SAMPLE(handles.quitfile)=[];         

    end

    filenamepd0=handles.filenamepd0;
    infile=handles.infile;
     
    filenamepd0(handles.quitfile)=[];
    infile(handles.quitfile)=[];
    handles.filenamepd0=filenamepd0;
    handles.infile=infile;
    handles.length=length(filenamepd0);
    guidata(hObject,handles)
end

% 
% --- Executes during object creation, after setting all properties.
function uitablesedimentology_CreateFcn(hObject, eventdata, handles)
% empty


% --- Executes when selected cell(s) is changed in inputfilename.
function inputfilename_CellSelectionCallback(hObject, eventdata, handles)
%define the number of file selected
handles.quitfile=eventdata.Indices(:,1);
guidata(hObject,handles)



% --- Executes when entered data in editable cell(s) in uitablesedimentology.
function uitablesedimentology_CellEditCallback(hObject, eventdata, handles)
% empty


% --- Executes when selected cell(s) is changed in uitablesedimentology.
function uitablesedimentology_CellSelectionCallback(hObject, eventdata, handles)
% empty     

% --------------------------------------------------------------------
function uitablesedimentology_ButtonDownFcn(hObject, eventdata, handles)
% empty
