function varargout = UI_SpFct(varargin)
% UI_SPFCT MATLAB code for UI_SpFct.fig
%      UI_SPFCT, by itself, creates a new UI_SPFCT or raises the existing
%      singleton*.
%
%      H = UI_SPFCT returns the handle to a new UI_SPFCT or the handle to
%      the existing singleton*.
%
%      UI_SPFCT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UI_SPFCT.M with the given input arguments.
%
%      UI_SPFCT('Property','Value',...) creates a new UI_SPFCT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UI_SpFct_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UI_SpFct_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UI_SpFct

% Last Modified by GUIDE v2.5 30-Oct-2018 14:49:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_SpFct_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_SpFct_OutputFcn, ...
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


% --- Executes just before UI_SpFct is made visible.
function UI_SpFct_OpeningFcn(hObject, eventdata, handles, varargin)

handles.amoy = 160;
handles.phi = 25;
handles.P = 0;
handles.alpha1 = [3.2/100 1];
handles.alpha0 = [0 1];
handles.c0 = 1.48;
handles.c1 = 0.64;
handles.rho0 = 1;
handles.rho1 = 1.85;
handles.NFmin = 0.5;
handles.NFmax = 4.5;
handles.NKmin = 0.5;
handles.NKmax = 4.5;
handles.Npts = 1000;
handles.Npks = 2;
handles.N = 10;
handles.filt = 0;
handles.seuil = 8;
handles.ceff = 0;
handles.SpFct = 0;
handles.f = 0;
handles.NormAxis = 0;
handles.fMin = handles.NFmin*handles.amoy/(4.*pi.*handles.amoy/1000);
handles.fMax = handles.NFmax*handles.amoy/(4.*pi.*handles.amoy/1000);
handles.load = 0;
handles.Search = 0;
handles.CalcD = 0;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = UI_SpFct_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PushReset.
function PushReset_Callback(hObject, eventdata, handles)
handles.amoy = 160;
handles.phi = 25;
handles.P = 0;
handles.alpha1 = [3.2/100 1];
handles.alpha0 = [0 1];
handles.NFmin = 0.5;
handles.NFmax = 4.5;
handles.NKmin = 0.5;
handles.NKmax = 4.5;
handles.Npts = 1000;
handles.Npks = 2;
handles.N = 10;
handles.filt = 0;
handles.seuil = 8;
handles.ceff = 0;
handles.SpFct = 0;
handles.f = 0;
handles.fMin = handles.NFmin*handles.amoy/(4.*pi.*handles.amoy/1000);
handles.fMax = handles.NFmax*handles.amoy/(4.*pi.*handles.amoy/1000);
handles.load = 0;
handles.Search = 0;
handles.CalcD = 0;

set(handles.EditNormFreqMin,'String',handles.NFmin)
set(handles.EditNormFreqMax,'String',handles.NFmax)
set(handles.EditNPeaks,'String',handles.Npks)
set(handles.EditN,'String',handles.N)
set(handles.EditNpts,'String',handles.Npts)
set(handles.EditNormKMin,'String',handles.NKmin)
set(handles.EditNormKMax,'String',handles.NKmax)
set(handles.EditAmoy,'String',handles.amoy)
set(handles.EditPhi,'String',handles.phi)
set(handles.EditAlpha0,'String',handles.alpha0(1).*100)
set(handles.EditAlpha1,'String',handles.alpha1(1).*100)
set(handles.EditPSD,'String',handles.P)
set(handles.EditfMin,'String',handles.EditfMin)
set(handles.EditfMax,'String',handles.EditfMax)
set(handles.Editc0,'String',handles.c0)
set(handles.Editc1,'String',handles.c1)
set(handles.Editrho0,'String',handles.rho0)
set(handles.Editrho1,'String',handles.rho1)


% --- Executes on button press in PushStart.
function PushStart_Callback(hObject, eventdata, handles)
[handles.SpFct,handles.ceff,handles.f] = calc_disp_mat(handles.amoy, ...
    handles.phi,handles.Npts,handles.NFmin ... 
    ,handles.NFmax,handles.N,handles.NKmin,handles.NKmax,handles.P ...
    ,handles.Npks,handles.alpha1,handles.alpha0,handles.filt, ...
    handles.seuil,handles.NormAxis,handles.Search,handles.CalcD, ...
    handles.c0,handles.c1,handles.rho0,handles.rho1);
guidata(hObject,handles)




function EditNormFreqMin_Callback(hObject, eventdata, handles)
handles.NFmin = str2double(get(hObject,'String'));
handles.fMin = handles.NFmin*handles.c0/(4.*pi.*handles.amoy/1000);
set(handles.EditfMin,'String',handles.fMin)
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function EditNormFreqMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNormFreqMax_Callback(hObject, eventdata, handles)
handles.NFmax = str2double(get(hObject,'String'));
handles.fMax = handles.NFmax*handles.c0/(4.*pi.*handles.amoy/1000);
set(handles.EditfMax,'String',handles.fMax)
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function EditNormFreqMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNPeaks_Callback(hObject, eventdata, handles)
handles.Npks = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during obje?ct creation, after setting all properties.
function EditNPeaks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditN_Callback(hObject, eventdata, handles)
handles.N = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function EditN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNpts_Callback(hObject, eventdata, handles)
handles.Npts = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditNpts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNormKMin_Callback(hObject, eventdata, handles)
handles.NKmin = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditNormKMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditNormKMax_Callback(hObject, eventdata, handles)
handles.NKmax = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditNormKMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditAmoy_Callback(hObject, eventdata, handles)
handles.amoy = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditAmoy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditPhi_Callback(hObject, eventdata, handles)
handles.phi = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditPhi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditAlpha0_Callback(hObject, eventdata, handles)
handles.alpha0 = [str2double(get(hObject,'String'))/100 1];
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditAlpha0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditAlpha1_Callback(hObject, eventdata, handles)
handles.alpha1 = [str2double(get(hObject,'String'))/100 1];
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditAlpha1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditPSD_Callback(hObject, eventdata, handles)
handles.P = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function EditPSD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSeuil_Callback(hObject, eventdata, handles)
handles.seuil = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function EditSeuil_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)

if handles.filt == 0
    
    handles.filt = 1;
    
else
    
    handles.filt = 0;
    
end

guidata(hObject,handles)


% --- Executes on button press in PushASCII.
function PushASCII_Callback(hObject, eventdata, handles)

    Chemin = uigetdir;
    
    dlmwrite([Chemin '\cEff.txt'],[handles.f.' handles.ceff'])
    dlmwrite([Chemin '\SpFct.txt'],handles.SpFct)


% --- Executes on button press in NormAxis.
function NormAxis_Callback(hObject, eventdata, handles)

if handles.NormAxis == 0
    
    handles.NormAxis = 1;
    
else
    
    handles.NormAxis = 0;
    
end

guidata(hObject,handles)



function EditfMin_Callback(hObject, eventdata, handles)
function EditfMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditfMax_Callback(hObject, eventdata, handles)
function EditfMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushLoad.
function PushLoad_Callback(hObject, eventdata, handles)

handles.load = 1;

[file,dir] = uigetfile('.mat');
load([dir file])

handles.amoy = amoy;
handles.phi = phi;
handles.P = P;
handles.alpha1 = alpha1;
handles.alpha0 = alpha0;


set(handles.EditAmoy,'String',handles.amoy)
set(handles.EditPhi,'String',handles.phi)
set(handles.EditAlpha0,'String',handles.alpha0(1))
set(handles.EditAlpha1,'String',handles.alpha1(1))
set(handles.EditPSD,'String',handles.P)


guidata(hObject,handles)


% --- Executes on button press in RadioSearch.
function RadioSearch_Callback(hObject, eventdata, handles)

if handles.Search == 0
    
    handles.Search = 1;
    
else
    
    handles.Search = 0;
    
end

guidata(hObject,handles)



function EditCalcD_Callback(hObject, eventdata, handles)

if handles.CalcD == 0
    
    handles.CalcD = 1;
    
else
    
    handles.CalcD = 0;
    
end

guidata(hObject,handles)


function Editc0_Callback(hObject, eventdata, handles)
handles.c0 = str2double(get(hObject,'String'));
guidata(hObject,handles)

function Editc0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Editc1_Callback(hObject, eventdata, handles)
handles.c1 = str2double(get(hObject,'String'));
guidata(hObject,handles)


function Editc1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Editrho0_Callback(hObject, eventdata, handles)
handles.rho0 = str2double(get(hObject,'String'));
guidata(hObject,handles)



function Editrho0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Editrho1_Callback(hObject, eventdata, handles)
handles.rho1 = str2double(get(hObject,'String'));
guidata(hObject,handles)


function Editrho1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
