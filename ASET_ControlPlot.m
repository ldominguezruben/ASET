function varargout = ASET_ControlPlot(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setting the cross section plot (Vel, Ms2, Gw and Gss). The
% range of contour plot, the colormap style, vertical exaggeration,
% smoothing axis and range axis.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_ControlPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_ControlPlot_OutputFcn, ...
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


% --- Executes just before ASET_ControlPlot is made visible.
function ASET_ControlPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for ASET_ControlPlot
handles.output = hObject;

guiparams = getappdata(0,'guiparams');

set(handles.figure1,'Name',['ASET: Plot Setting '], ...
'DockControls','off');

handles.first=1;
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ASET_ControlPlot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in applypush.
function applypush_Callback(hObject, eventdata, handles,varargin)

cof.crosssel=get(handles.varselect,'Value');

if handles.first==1
    Method=1;%Default graph

    handles.first=0;%change option
    cof.verticalsmooth=0;
    cof.horizontalsmooth=0;
    [handles.sal]=ASET_PlotXS(cof,Method);
    % Update handles structure
    guidata(hObject, handles);

    load_init(handles)

else
    Method=2;%modify

    %load variables init
    cof.minxaxis=str2double(get(handles.minxaxis,'String'));
    cof.maxxaxis=str2double(get(handles.maxxaxis,'String'));
    cof.minyaxis=str2double(get(handles.minyaxis,'String'));
    cof.maxyaxis=str2double(get(handles.maxyaxis,'String'));
    cof.minvel=str2double(get(handles.minvel,'String'));
    cof.maxvel=str2double(get(handles.maxvel,'String'));

    cof.verticalexaggeration=str2double(get(handles.verticalexaggeration,'String'));
    cof.verticalsmooth=str2double(get(handles.verticalsmooth,'String'));
    cof.horizontalsmooth=str2double(get(handles.horizontalsmooth,'String'));

    % Colormap choices

    idx=get(handles.colormappopup,'Value');

    switch idx
        case 1
            cof.colormaps_mcs   = 'Jet';
        case 2
            cof.colormaps_mcs   = 'HSV';
        case 3
            cof.colormaps_mcs   = 'Hot';
        case 4
            cof.colormaps_mcs   = 'Cool';
        case 5
            cof.colormaps_mcs   = 'Spring';
        case 6
            cof.colormaps_mcs   = 'Summer';
        case 7
            cof.colormaps_mcs   = 'Autumn';    
        case 8
            cof.colormaps_mcs   = 'Winter';
        case 9
            cof.colormaps_mcs   = 'Gray';
        case 10
            cof.colormaps_mcs   = 'Bone';
        case 11
            cof.colormaps_mcs   = 'Copper';
        case 12
            cof.colormaps_mcs   = 'Pink';
    end

    [handles.sal]=ASET_PlotXS(cof,Method);
    % Update handles structure
    guidata(hObject, handles);
end


function minvel_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function minvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxvel_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function maxvel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxvel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minxaxis_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function minxaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minxaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxxaxis_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function maxxaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxxaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minyaxis_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function minyaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minyaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxyaxis_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function maxyaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxyaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colormappopup.
function colormappopup_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function colormappopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colormappopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function verticalexaggeration_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function verticalexaggeration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verticalexaggeration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function verticalsmooth_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function verticalsmooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verticalsmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function horizontalsmooth_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function horizontalsmooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to horizontalsmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_init(handles)

%load variables init
set(handles.minxaxis,'String',num2str(round(handles.sal.xlim(1)*1000)/1000))
set(handles.maxxaxis,'String',num2str(round(handles.sal.xlim(2)*100)/100))
set(handles.minyaxis,'String',num2str(round(handles.sal.ylim(1)*1000)/1000))
set(handles.maxyaxis,'String',num2str(round(handles.sal.ylim(2)*100)/100))
set(handles.minvel,'String',num2str(handles.sal.zmin))
set(handles.maxvel,'String',num2str(handles.sal.zmax))
set(handles.maxvel,'String',num2str(handles.sal.zmax))

set(handles.verticalexaggeration,'String',num2str(handles.sal.exaggeration))
set(handles.verticalsmooth,'String',num2str(handles.sal.smoothvertical))
set(handles.horizontalsmooth,'String',num2str(handles.sal.smoothhorizontal))


% --- Executes on selection change in varselect.
function varselect_Callback(hObject, eventdata, handles)
handles.first=1;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function varselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over varselect.
function varselect_ButtonDownFcn(hObject, eventdata, handles)
% empty


% --- Executes on key press with focus on varselect and none of its controls.
function varselect_KeyPressFcn(hObject, eventdata, handles)
% empty
