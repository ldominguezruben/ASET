function varargout = ASET_ADCPData(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function acquires the data from the .mmt file that you read and plot 
% them allowing you to modify and re-edit them. You also read the draft and
% ADCP frequency

%Dominguez Ruben 5/10/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_ADCPData_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_ADCPData_OutputFcn, ...
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


% --- Executes just before ASET_ADCPData is made visible.
function ASET_ADCPData_OpeningFcn(hObject, eventdata, handles, varargin)


set(handles.figure1,'Name',['ASET: ADCP and Measurement Settings '], ...
    'DockControls','off')

% Input variables for rio Grande 1200kHz (Teledyne)
if varargin{1}.Calibration==1 %Calibration Module 
    if varargin{1}.first==0; %firs read
        %Set data for Calibration
        handles.noenable=1;
        handles.noenablekcEr=handles.noenable;
        handles.first=varargin{1}.first;
        handles.R.multi=varargin{1}.R.multi;
        handles.length=varargin{1}.length;
        handles.Calibration=varargin{1}.Calibration;
        set(handles.kcEr,'Enable','off');
        set(handles.readcalibrationfile,'Enable','off');
        guidata(hObject,handles)
    elseif varargin{1}.first==1; %second read
        %read and write data
        set(handles.rssi_beam1,'String',num2str(varargin{1}.rssi_beam1(1)));
        set(handles.rssi_beam2,'String',num2str(varargin{1}.rssi_beam2(1)));
        set(handles.rssi_beam3,'String',num2str(varargin{1}.rssi_beam3(1)));
        set(handles.rssi_beam4,'String',num2str(varargin{1}.rssi_beam4(1)));
        set(handles.frequency, 'String', varargin{1}.Inst.freq(1))
        set(handles.adcptype, 'Value', varargin{1}.adcp)
        set(handles.ccoef,'String',num2str(varargin{1}.ccoef));
        set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))
        %Set data for Calibration
        handles.noenable=1;
        handles.noenablekcEr=handles.noenable;
        handles.first=varargin{1}.first;
        handles.R.multi=varargin{1}.R.multi;
        handles.length=varargin{1}.length;
        handles.Calibration=varargin{1}.Calibration;
        set(handles.kcEr,'Enable','off');
        set(handles.readcalibrationfile,'Enable','off');
        guidata(hObject,handles)
    end
    
    
elseif varargin{1}.Calibration==0 % simple or multiple select file
    
    handles.noenable=0;

    if varargin{1}.R.multi==0
        set(handles.draft,'String',num2str(varargin{1}.R.draft))
        P.Ac.FREQ=varargin{1}.Inst.freq(1);
        set(handles.frequency,'String',num2str(P.Ac.FREQ));
    elseif varargin{1}.R.multi==1
        set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))    
        P.Ac.FREQ=varargin{1}.Inst.freq(1);
        set(handles.frequency,'String',num2str(P.Ac.FREQ));
    end

    
     if P.Ac.FREQ==600
        P.Ac.beam_1=0.4126;%beam coefficient 600KHz FICH UNL
        P.Ac.beam_2=0.4027;%beam coefficient 600KHz FICH UNL
        P.Ac.beam_3=0.4;%beam coefficient 600KHz FICH UNL
        P.Ac.beam_4=0.4028;%beam coefficient 600KHz FICH UNL
        num=1;
        ccoef=-139.09;
    elseif P.Ac.FREQ==1200
        P.Ac.beam_1=0.3909;%beam coefficient 1200KHz FICH UNL
        P.Ac.beam_2=0.4094;%beam coefficient 1200KHz FICH UNL
        P.Ac.beam_3=0.4061;%beam coefficient 1200KHz FICH UNL
        P.Ac.beam_4=0.412;%beam coefficient 1200KHz FICH UNL
        num=2;
        ccoef=-129.09;
    end

    if varargin{1}.first==0; %primera lectura 
        if varargin{1}.R.multi==0%single select
            %read default data
            set(handles.rssi_beam1,'String',num2str(P.Ac.beam_1));
            set(handles.rssi_beam2,'String',num2str(P.Ac.beam_2));
            set(handles.rssi_beam3,'String',num2str(P.Ac.beam_3));
            set(handles.rssi_beam4,'String',num2str(P.Ac.beam_4));
            set(handles.adcptype, 'Value', num)
            set(handles.ccoef,'String',num2str(ccoef));
            set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))
        elseif varargin{1}.R.multi==1;%multiselect
            set(handles.rssi_beam1,'String',num2str(P.Ac.beam_1));
            set(handles.rssi_beam2,'String',num2str(P.Ac.beam_2));
            set(handles.rssi_beam3,'String',num2str(P.Ac.beam_3));
            set(handles.rssi_beam4,'String',num2str(P.Ac.beam_4));
            set(handles.adcptype, 'Value', num)
            set(handles.ccoef,'String',num2str(ccoef))
            set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))
        end
    elseif varargin{1}.first==1; %seconda read
        %Read previous data
        if varargin{1}.R.multi==0;%single select
            set(handles.rssi_beam1,'String',num2str(varargin{1}.rssi_beam1));
            set(handles.rssi_beam2,'String',num2str(varargin{1}.rssi_beam2));
            set(handles.rssi_beam3,'String',num2str(varargin{1}.rssi_beam3));
            set(handles.rssi_beam4,'String',num2str(varargin{1}.rssi_beam4));
            set(handles.adcptype, 'Value', num)
            set(handles.ccoef,'String',num2str(ccoef));
            set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))
        elseif varargin{1}.R.multi==1;%multi select
            set(handles.rssi_beam1,'String',num2str(varargin{1}.rssi_beam1(1)));
            set(handles.rssi_beam2,'String',num2str(varargin{1}.rssi_beam2(1)));
            set(handles.rssi_beam3,'String',num2str(varargin{1}.rssi_beam3(1)));
            set(handles.rssi_beam4,'String',num2str(varargin{1}.rssi_beam4(1)));
            set(handles.adcptype, 'Value', num)
            set(handles.ccoef,'String',num2str(ccoef));
            set(handles.draft,'String',num2str(varargin{1}.R.draft(1)))
        end
    end

    if varargin{1}.first==1;%second read
        set(handles.kcEr,'String',num2str(varargin{1}.kcEr));
    end
        set(handles.kcEr,'Enable','on');

    % Save the data 
    handles.R.multi=varargin{1}.R.multi;
    handles.first=varargin{1}.first;
    handles.Freq=P.Ac.FREQ;
    handles.length=varargin{1}.length;
    handles.Calibration=varargin{1}.Calibration;
    guidata(hObject,handles)
        
end
   % Choose default command line output for ASET_ADCPData
    handles.output = hObject;
    guidata(hObject,handles)
        

% --- Outputs from this function are returned to the command line.
function varargout = ASET_ADCPData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%
%ADCP Data
%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function adcptype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adcptype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in adcptype.
function adcptype_Callback(hObject, eventdata, handles)

freq=get(handles.adcptype,'Value');

 if freq==1%rio Grande 600kHz
    P.Ac.beam_1=0.4126;%beam coefficient 600KHz
    P.Ac.beam_2=0.4027;%beam coefficient 600KHz
    P.Ac.beam_3=0.4;%beam coefficient 600KHz
    P.Ac.beam_4=0.4028;%beam coefficient 600KHz
    frequenc=600;
    ccoef=-139.09;
elseif freq==2%rio Grande 1200kHz
    P.Ac.beam_1=0.3909;%beam coefficient 1200KHz
    P.Ac.beam_2=0.4094;%beam coefficient 1200KHz
    P.Ac.beam_3=0.4061;%beam coefficient 1200KHz
    P.Ac.beam_4=0.412;%beam coefficient 1200KHz
    frequenc=1200;
    ccoef=-129.09;
elseif freq==3%RioPro 1200kHz
    P.Ac.beam_1=[];
    P.Ac.beam_2=[];
    P.Ac.beam_3=[];
    P.Ac.beam_4=[];
    frequenc=1200;
    ccoef=-131.36;
 elseif freq==4%RiverPro 1200kHz
     P.Ac.beam_1=[];
     P.Ac.beam_2=[];
     P.Ac.beam_3=[];
     P.Ac.beam_4=[];
     frequenc=1200;
     ccoef=-128.09;
 elseif freq==5%RiverRay 600kHz
     P.Ac.beam_1=0.6;
     P.Ac.beam_2=0.6;
     P.Ac.beam_3=0.6;
     P.Ac.beam_4=0.6;
     frequenc=600;
     ccoef=-138.02;
end
                           
%Read default data
set(handles.rssi_beam1,'String',num2str(P.Ac.beam_1));
set(handles.rssi_beam2,'String',num2str(P.Ac.beam_2));
set(handles.rssi_beam3,'String',num2str(P.Ac.beam_3));
set(handles.rssi_beam4,'String',num2str(P.Ac.beam_4));
set(handles.ccoef,'String',num2str(ccoef));
set(handles.frequency,'String',frequenc);

function rssi_beam1_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function rssi_beam1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rssi_beam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rssi_beam2_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function rssi_beam2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rssi_beam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rssi_beam3_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function rssi_beam3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rssi_beam3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rssi_beam4_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function rssi_beam4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rssi_beam4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function draft_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function draft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to draft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function kcEr_Callback(hObject, eventdata, handles)
% Empty


% --- Executes during object creation, after setting all properties.
function kcEr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kcEr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calibration Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in readcalibrationfile.
function readcalibrationfile_Callback(hObject, eventdata, handles)
[filename,pathname] = uigetfile({'*.mat';'*.*'},'Select .mat the file');
  
if ischar(filename) 
    %Read .mat file (from VMT USGS)
    vars=load(fullfile(pathname,filename));
else
end

handles.Cal=vars.Cal;
set(handles.kcEr,'String',num2str(handles.Cal.kcEr_mean));
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in saveadcpdata.
function saveadcpdata_Callback(hObject, eventdata, handles)

if handles.noenable==1
    if isempty(str2num(get(handles.rssi_beam1,'String'))) 
        wrong_data=warndlg('Please complete the Beam 1 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam2,'String')))
        wrong_data=warndlg('Please complete the Beam 2 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam3,'String')))
        wrong_data=warndlg('Please complete the Beam 3 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam4,'String')))
        wrong_data=warndlg('Please complete the Beam 4 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.draft,'String')))
        wrong_data=warndlg('Please complete the Draft box',...
        'WARNING'); 
        h=0;
    else
        h=1;
    end
       
else
    if isempty(str2num(get(handles.rssi_beam1,'String'))) 
        wrong_data=warndlg('Please complete the Beam 1 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam2,'String')))
        wrong_data=warndlg('Please complete the Beam 2 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam3,'String')))
        wrong_data=warndlg('Please complete the Beam 3 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.rssi_beam4,'String')))
        wrong_data=warndlg('Please complete the Beam 4 box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.draft,'String')))
        wrong_data=warndlg('Please complete the Draft box',...
        'WARNING'); 
        h=0;
    elseif isempty(str2num(get(handles.kcEr,'String')))
        wrong_data=warndlg('Please complete the kcEr box',...
        'WARNING'); 
        h=0;
    else
        h=1;
    end
end
    
% Define the frequency of ADCP equipment
freq=get(handles.adcptype,'Value');

 if freq==1
    varargout.Inst.freq(1)=600; 
 elseif freq==2
    varargout.Inst.freq(1)=1200;
 end
 
 % ADCP type
 varargout.typeADCP=handles.adcptype.String(handles.adcptype.Value);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
if handles.first==0%first read
    for i=1:handles.length%single or multiple bean
        rssi_beam1=get(handles.rssi_beam1,'String');
        varargout.rssi_beam1(i)=str2num(rssi_beam1);

        rssi_beam2=get(handles.rssi_beam2,'String');
        varargout.rssi_beam2(i)=str2num(rssi_beam2);

        rssi_beam3=get(handles.rssi_beam3,'String');
        varargout.rssi_beam3(i)=str2num(rssi_beam3);

        rssi_beam4=get(handles.rssi_beam4,'String');
        varargout.rssi_beam4(i)=str2num(rssi_beam4);

        draft=get(handles.draft,'String');
        varargout.R.draft(i)=str2num(draft);
    end
elseif handles.first==1%second read
    varargout = getappdata(0,'ADCPdata');%read varargout

    if handles.length==1%single select
        rssi_beam1=get(handles.rssi_beam1,'String');
        varargout.rssi_beam1=str2num(rssi_beam1);

        rssi_beam2=get(handles.rssi_beam2,'String');
        varargout.rssi_beam2=str2num(rssi_beam2);

        rssi_beam3=get(handles.rssi_beam3,'String');
        varargout.rssi_beam3=str2num(rssi_beam3);

        rssi_beam4=get(handles.rssi_beam4,'String');
        varargout.rssi_beam4=str2num(rssi_beam4);

        draft=get(handles.draft,'String');
        varargout.R.draft=str2num(draft);

    elseif handles.length>1
        for i=1:handles.length% multiple select
            rssi_beam1=get(handles.rssi_beam1,'String');
            varargout.rssi_beam1(i)=str2num(rssi_beam1);

            rssi_beam2=get(handles.rssi_beam2,'String');
            varargout.rssi_beam2(i)=str2num(rssi_beam2);

            rssi_beam3=get(handles.rssi_beam3,'String');
            varargout.rssi_beam3(i)=str2num(rssi_beam3);

            rssi_beam4=get(handles.rssi_beam4,'String');
            varargout.rssi_beam4(i)=str2num(rssi_beam4);

            draft=get(handles.draft,'String');
            varargout.R.draft(i)=str2num(draft);
        end
    end
end
        
%Save data
kcEr=get(handles.kcEr,'String');
varargout.kcEr=str2num(kcEr);
handles.first=1;
varargout.adcp=handles.adcptype.Value;
varargout.R.multi=handles.R.multi;
varargout.first=handles.first;
varargout.length=handles.length;
varargout.Calibration=handles.Calibration;
varargout.ccoef=get(handles.ccoef,'String');
varargout.saveADCPData=1;
setappdata(0,'ADCPdata',varargout)

                    
%Close the GUI if is complete the data
 if h==0
 elseif h==1
    close(handles.figure1)
 end
 

% %%%%%%%%%%%%%%%%%%%%%%%
% %Extra Function
% %%%%%%%%%%%%%%%%%%%%%%%
function frequency_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ccoef_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function ccoef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ccoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
