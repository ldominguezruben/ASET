function varargout = ASET_VerticalProfilePlot(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plot the vertical profile of velocity and concentration.

% by Dominguez Ruben, L. FICH-UNL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASET_VerticalProfilePlot_OpeningFcn, ...
                   'gui_OutputFcn',  @ASET_VerticalProfilePlot_OutputFcn, ...
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


% --- Executes just before ASET_VerticalProfilePlot is made visible.
function ASET_VerticalProfilePlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ASET_VerticalProfilePlot (see VARARGIN)

% Choose default command line output for ASET_VerticalProfilePlot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'Name',['ASET: Vertical Profile'],...
    'DockControls','off')

% UIWAIT makes ASET_VerticalProfilePlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ASET_VerticalProfilePlot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plotverticalpro.
function plotverticalpro_Callback(hObject, eventdata, handles)

vertical=get(handles.verticalprofilelist,'Value');

%Load data
guiparams = getappdata(0,'guiparams');
V   = guiparams.V;
Array   = guiparams.Array;

switch vertical
    case 1
        figure
        plot(Array.Css*1000,Array.Depth_Css,'ok');
        set(gca,'ydir','reverse')
        ylim([0 max(max(V.mcsBed))])   
        title('ASET: Ms2 Vertical Profile [mg/l]','Fontsize',18)
        ylabel('Depth [m]','Fontsize',16);xlabel('Ms2 [mg/l]','Fontsize',16)
    case 2
        figure
        plot(Array.uMag,Array.Depth_uMag,'ok');
        set(gca,'ydir','reverse')
        ylim([0 max(max(V.mcsBed))])
        title('ASET: Velocity Vertical Profile [cm/s]','Fontsize',18)
        ylabel('Depth [m]','Fontsize',16);xlabel('Velocity Magnitude [cm/s]','Fontsize',16)
end

% --- Executes on selection change in verticalprofilelist.
function verticalprofilelist_Callback(hObject, eventdata, handles)
% empty


% --- Executes during object creation, after setting all properties.
function verticalprofilelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verticalprofilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
