function varargout = SpatialSummationConfig(varargin)
% SPATIALSUMMATIONCONFIG M-file for SpatialSummationConfig.fig
%      SPATIALSUMMATIONCONFIG, by itself, creates a new SPATIALSUMMATIONCONFIG or raises the existing
%      singleton*.
%
%      H = SPATIALSUMMATIONCONFIG returns the handle to a new SPATIALSUMMATIONCONFIG or the handle to
%      the existing singleton*.
%
%      SPATIALSUMMATIONCONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPATIALSUMMATIONCONFIG.M with the given input arguments.
%
%      SPATIALSUMMATIONCONFIG('Property','Value',...) creates a new SPATIALSUMMATIONCONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpatialSummationConfig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpatialSummationConfig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpatialSummationConfig

% Last Modified by GUIDE v2.5 29-Aug-2016 10:35:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SpatialSummationConfig_OpeningFcn, ...
    'gui_OutputFcn',  @SpatialSummationConfig_OutputFcn, ...
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


% --- Executes just before SpatialSummationConfig is made visible.
function SpatialSummationConfig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpatialSummationConfig (see VARARGIN)

% Choose default command line output for SpatialSummationConfig
global StimParams VideoParams;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpatialSummationConfig wait for user response (see UIRESUME)
% uiwait(handles.main_CFG_figure);
hAomControl = getappdata(0,'hAomControl');
exp = getappdata(hAomControl, 'exp');

if exist('lastBackgroundExchangeCFG.mat','file')==2
    loadLastValues(handles,1);
else
    loadLastValues(handles,0);
end

set(handles.stimpath, 'String', [cd '\tempStimulus\']);
set(handles.ok_button, 'Enable', 'on');
set(handles.gain_lbl, 'Visible', 'On');
set(handles.gain, 'Visible', 'On');


% --- Outputs from this function are returned to the command line.
function varargout = SpatialSummationConfig_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function initials_Callback(hObject, eventdata, handles)
% hObject    handle to initials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initials as text
%        str2double(get(hObject,'String')) returns contents of initials as a double
user_entry = get(hObject,'string');
if get(handles.auto_prefix, 'Value') == 1
    set(handles.vidprefix, 'String', user_entry)
elseif get(handles.auto_prefix, 'Value') == 0
    %do nothing
end


% --- Executes during object creation, after setting all properties.
function initials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function pupilsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupilsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function fieldsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function videodur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to videodur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function vidprefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vidprefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function presentdur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to presentdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function npresent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npresent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function thresholdGuess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function priorSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to priorSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_stimpath.
function set_stimpath_Callback(hObject, eventdata, handles)
% hObject    handle to set_stimpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stimpath = uigetdir('path','Select directory containing the stimuli');
set(handles.stimpath, 'String', [stimpath '\']);
if stimpath == 0;
    %do nothing
else
    set(handles.ok_button, 'Enable', 'on');
    hAomControl = getappdata(0,'hAomControl');
    setappdata(hAomControl, 'stimpath', stimpath);
end


% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SYSPARAMS StimParams VideoParams ExpCfgParams;

hAomControl = getappdata(0,'hAomControl');
stimpath = getappdata(hAomControl, 'stimpath');
CFG.ok = 1;
CFG.initials = get(handles.initials, 'String');
CFG.pupilsize = get(handles.pupilsize, 'String');
CFG.fieldsize = str2double(get(handles.fieldsize, 'String'));
CFG.presentdur = str2double(get(handles.presentdur, 'String'));
CFG.videodur = str2double(get(handles.videodur, 'String'));
CFG.vidprefix = get(handles.vidprefix, 'String');
CFG.kb_next = get(handles.kb_next, 'String');
CFG.kb_yes = get(handles.kb_yes, 'String');
CFG.kb_no = get(handles.kb_no, 'String');
CFG.kb_repeat = get(handles.kb_repeat, 'String');
CFG.record = 1;

if get(handles.quest_radio, 'Value') == 1;
    CFG.method = 'q';
    CFG.quest_radio = 1; 
end

if get(handles.yes_no_radio, 'Value') == 1;
    CFG.subject_response = 'y';
    CFG.yes_no_radio = 1;
    CFG.two_afc_radio = 0;
elseif get(handles.two_afc_radio, 'Value') == 1;
    CFG.subject_response = '2';
    CFG.yes_no_radio = 0;
    CFG.two_afc_radio = 1;
else
end

if get(handles.circle_button, 'Value') == 1;
    CFG.stim_shape = 'Circle';
    CFG.circle_button = 1;
    CFG.square_button = 0;
elseif get(handles.square_button, 'Value') == 1;
    CFG.stim_shape = 'Square';
    CFG.circle_button = 0;
    CFG.square_button = 1;
else
end

CFG.npresent = str2double(get(handles.npresent, 'String'));
CFG.gain = str2double(get(handles.gain, 'String'));
CFG.beta = str2double(get(handles.beta, 'String'));
CFG.delta = str2double(get(handles.delta, 'String'));
CFG.pCorrect = str2double(get(handles.pCorrect, 'String'));
CFG.thresholdGuess = str2double(get(handles.thresholdGuess, 'String'));
CFG.priorSD = str2double(get(handles.priorSD, 'String'));
CFG.stimpath = get(handles.stimpath,'String');
CFG.green_x_offset = str2double(get(handles.green_x_offset,'String'));
CFG.green_y_offset = str2double(get(handles.green_y_offset,'String'));
CFG.red_x_offset = str2double(get(handles.red_x_offset,'String'));
CFG.red_y_offset = str2double(get(handles.red_y_offset,'String'));
CFG.red_stim_color = get(handles.red_stim_radio, 'Value');
CFG.green_stim_color = get(handles.green_stim_radio, 'Value');
CFG.stim_midpoint = str2double(get(handles.stim_midpoint, 'String'));
CFG.num_stims_riccos = str2double(get(handles.num_stims_riccos, 'String'));
CFG.frame_rate = str2double(get(handles.frame_rate, 'String'));
CFG.presentdur_frames = str2double(get(handles.presentdur_frames, 'String'));
CFG.interleave_check = get(handles.interleave_check, 'Value');
CFG.low_pCorrect = str2double(get(handles.low_pCorrect, 'String'));
CFG.high_pCorrect = str2double(get(handles.high_pCorrect, 'String'));
CFG.gamepad_radio = get(handles.gamepad_radio, 'Value');
CFG.keyboard_radio = get(handles.keyboard_radio, 'Value');

% CFG.last_offset_filename = ExpCfgParams.last_offset_filename;


ExpCfgParams.initials = CFG.initials;
ExpCfgParams.pupilsize = CFG.pupilsize;
ExpCfgParams.presentdur = CFG.presentdur;
ExpCfgParams.fieldsize = CFG.fieldsize; % Rastersize in degrees
ExpCfgParams.videodur = CFG.videodur; % Video recording time in seconds
ExpCfgParams.quest_radio = CFG.quest_radio;
ExpCfgParams.priorSD = CFG.priorSD;
ExpCfgParams.text17 = 1;
ExpCfgParams.text19 = 1;
ExpCfgParams.pCorrect = CFG.pCorrect;
ExpCfgParams.text15 = 1;
ExpCfgParams.thresholdGuess = CFG.thresholdGuess;
ExpCfgParams.text14 = 1;
ExpCfgParams.beta = CFG.beta;
ExpCfgParams.text16 = 1;
ExpCfgParams.delta = CFG.delta;
ExpCfgParams.circle_button = CFG.circle_button;
ExpCfgParams.square_button = CFG.square_button;
ExpCfgParams.yes_no_radio = CFG.yes_no_radio;
ExpCfgParams.two_afc_radio = CFG.two_afc_radio;
ExpCfgParams.npresent = CFG.npresent;
ExpCfgParams.gain = CFG.gain;
ExpCfgParams.green_x_offset = CFG.green_x_offset;
ExpCfgParams.green_y_offset = CFG.green_y_offset;
ExpCfgParams.red_x_offset = CFG.red_x_offset;
ExpCfgParams.red_y_offset = CFG.red_y_offset;
ExpCfgParams.red_stim_color = CFG.red_stim_color;
ExpCfgParams.green_stim_color = CFG.green_stim_color;
ExpCfgParams.kb_yes = CFG.kb_yes;
ExpCfgParams.kb_no = CFG.kb_no;
ExpCfgParams.kb_repeat = CFG.kb_repeat;
ExpCfgParams.kb_next = CFG.kb_next;
ExpCfgParams.stim_midpoint = CFG.stim_midpoint;
ExpCfgParams.num_stims_riccos = CFG.num_stims_riccos;
ExpCfgParams.frame_rate = CFG.frame_rate;
ExpCfgParams.presentdur_frames = CFG.presentdur_frames;
ExpCfgParams.interleave_check = CFG.interleave_check;
ExpCfgParams.low_pCorrect = CFG.low_pCorrect;
ExpCfgParams.high_pCorrect = CFG.high_pCorrect;
ExpCfgParams.gamepad_radio = CFG.gamepad_radio;
ExpCfgParams.keyboard_radio = CFG.keyboard_radio;


setappdata(hAomControl, 'CFG', CFG);
save('lastBackgroundExchangeCFG.mat', 'ExpCfgParams','CFG')
% save(fullfile('Experiments', 'lastBackgroundExchangeCFG.mat'), 'ExpCfgParams', 'CFG'); 
close;

% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hAomControl = getappdata(0,'hAomControl');
CFG.ok = 0;
setappdata(hAomControl, 'CFG', CFG);
close;


% --- Executes during object creation, after setting all properties.
function pCorrect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in method_radio_panel.
function method_radio_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in method_radio_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global SYSPARAMS StimParams VideoParams ExpCfgParams;

if get(handles.quest_radio, 'Value') == 1;
    set(handles.priorSD, 'Visible', 'on');
    set(handles.priorSD, 'String', 1);
    set(handles.text14, 'Visible', 'on', 'String', 'Beta');
    set(handles.beta, 'Visible', 'on');
    set(handles.text16, 'Visible', 'on');
    set(handles.delta, 'Visible', 'on');
    set(handles.text19, 'Visible', 'on');
    set(handles.pCorrect, 'Visible', 'on');
    set(handles.text17, 'String', 'Prior SD', 'Visible', 'on');
    set(handles.text13, 'Visible', 'on');
    set(handles.npresent, 'Visible', 'on', 'Enable', 'on');
    set(handles.quest_params, 'Title', 'QUEST Parameters');
    set(handles.text14, 'Visible', 'on', 'String', 'Beta');
    set(handles.beta, 'Visible', 'on', 'String', 3.5);
    set(handles.text16, 'Visible', 'on', 'String', 'Delta');
%     set(handles.delta, 'Visible', 'on', 'String', 0.01);
%     set(handles.repeat_string, 'Visible', 'off');
%     set(handles.numrepeat, 'Visible', 'off');
% elseif get(handles.fourtwo_radio, 'Value') == 1;
%     set(handles.priorSD, 'Visible', 'off');
%     set(handles.text14, 'Visible', 'off');
%     set(handles.beta, 'Visible', 'off');
%     set(handles.text16, 'Visible', 'off');
%     set(handles.delta, 'Visible', 'off');
%     set(handles.text19, 'Visible', 'off');
%     set(handles.pCorrect, 'Visible', 'off');
%     set(handles.text17, 'Visible', 'off');
%     set(handles.text13, 'Visible', 'on');
%     set(handles.npresent, 'Value', 40,'Enable', 'off');
%     set(handles.quest_params, 'Title', 'Staircase Parameters');
%     set(handles.npresent, 'Value', 40);
%     set(handles.repeat_string, 'Visible', 'on');
%     set(handles.numrepeat, 'Visible', 'on');
end


% --- Executes when selected object is changed in method_uipanel.
function method_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in method_uipanel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.yes_no_radio,'Value')==1;
    set(handles.pCorrect, 'String', 50);
elseif get(handles.two_afc_radio,'Value')==1;
    set(handles.pCorrect, 'String', 75);
else
    set(handles.pCorrect,'String', 62.5);
end


% --- Executes during object creation, after setting all properties.
function stim_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function green_y_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_y_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function green_x_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_x_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function red_y_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_y_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function red_x_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_x_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function loadLastValues (handles,last)

global SYSPARAMS StimParams VideoParams ExpCfgParams;

if last==1
    load('lastBackgroundExchangeCFG.mat');
else %use defaults
    ExpCfgParams.initials = 'Sample';
    ExpCfgParams.pupilsize = 6.0;
    ExpCfgParams.presentdur = 200;
    ExpCfgParams.fieldsize = 420; % Rastersize in ppd
    ExpCfgParams.videodur = 1; % Video recording time in seconds
    ExpCfgParams.quest_radio = 1;
    ExpCfgParams.priorSD = 1;
    ExpCfgParams.text17 = 1;
    ExpCfgParams.text19 = 1;
    ExpCfgParams.pCorrect = 78;
    ExpCfgParams.text15 = 1;
    ExpCfgParams.thresholdGuess = 0.5;
    ExpCfgParams.text14 = 1;
    ExpCfgParams.beta = 3.5;
    ExpCfgParams.text16 = 1;
    ExpCfgParams.delta = 0.01;
    ExpCfgParams.circle_button = 1;
    ExpCfgParams.square_button = 0;
    ExpCfgParams.yes_no_radio = 1;
    ExpCfgParams.two_afc_radio = 0;
    ExpCfgParams.npresent = 40;
    ExpCfgParams.gain = 1;
    ExpCfgParams.green_x_offset = 0;
    ExpCfgParams.green_y_offset = 0;
    ExpCfgParams.red_x_offset = 0;
    ExpCfgParams.red_y_offset = 0;
    ExpCfgParams.red_stim_color = 1;
    ExpCfgParams.green_stim_color = 0;
    ExpCfgParams.kb_next = 'space';
    ExpCfgParams.kb_yes = 'rightarrow';
    ExpCfgParams.kb_no = 'leftarrow';
    ExpCfgParams.kb_repeat = 'uparrow';
    ExpCfgParams.stim_midpoint = 16;
    ExpCfgParams.num_stims_riccos = 7;
    ExpCfgParams.frame_rate = 30;
    ExpCfgParams.presentdur_frames = round(ExpCfgParams.frame_rate*(ExpCfgParams.presentdur/1000));
    ExpCfgParams.gamepad_radio = 0;
    ExpCfgParams.keyboard_radio = 1;
    ExpCfgParams.interleave_check = 0;
    ExpCfgParams.low_pCorrect = 50;
    ExpCfgParams.high_pCorrect = 90;
    ExpCfgParams.last_offset_filename = 'test';
    
end

set(handles.initials, 'String', ExpCfgParams.initials);
set(handles.pupilsize, 'String', ExpCfgParams.pupilsize);
set(handles.fieldsize, 'String', ExpCfgParams.fieldsize);
set(handles.presentdur, 'String', ExpCfgParams.presentdur);
set(handles.videodur, 'String', ExpCfgParams.videodur);
set(handles.quest_radio, 'Value', ExpCfgParams.quest_radio);
set(handles.priorSD, 'String', ExpCfgParams.priorSD);
set(handles.pCorrect, 'String', ExpCfgParams.pCorrect);
set(handles.thresholdGuess, 'String', ExpCfgParams.thresholdGuess);
set(handles.beta, 'String', ExpCfgParams.beta);
set(handles.delta, 'String', ExpCfgParams.delta);
set(handles.circle_button, 'Value', ExpCfgParams.circle_button);
set(handles.square_button, 'Value', ExpCfgParams.square_button);
set(handles.yes_no_radio, 'Value', ExpCfgParams.yes_no_radio);
set(handles.two_afc_radio, 'Value', ExpCfgParams.two_afc_radio);
set(handles.npresent, 'String', ExpCfgParams.npresent);
set(handles.gain, 'String', ExpCfgParams.gain);
set(handles.green_x_offset, 'String', ExpCfgParams.green_x_offset);
set(handles.green_y_offset, 'String', ExpCfgParams.green_y_offset);
set(handles.red_x_offset, 'String', ExpCfgParams.red_x_offset);
set(handles.red_y_offset, 'String', ExpCfgParams.red_y_offset);
set(handles.red_stim_radio, 'Value', ExpCfgParams.red_stim_color);
set(handles.green_stim_radio, 'Value', ExpCfgParams.green_stim_color);
set(handles.kb_next, 'String', ExpCfgParams.kb_next);
set(handles.kb_yes, 'String', ExpCfgParams.kb_yes);
set(handles.kb_no, 'String', ExpCfgParams.kb_no);
set(handles.kb_repeat, 'String', ExpCfgParams.kb_repeat);
set(handles.stim_midpoint, 'String', ExpCfgParams.stim_midpoint);
set(handles.num_stims_riccos, 'String', ExpCfgParams.num_stims_riccos);
set(handles.frame_rate, 'String', ExpCfgParams.frame_rate);
set(handles.presentdur_frames, 'String', ExpCfgParams.presentdur_frames);
set(handles.interleave_check, 'Value', ExpCfgParams.interleave_check);
set(handles.low_pCorrect, 'String', ExpCfgParams.low_pCorrect);
set(handles.high_pCorrect, 'String', ExpCfgParams.high_pCorrect);
set(handles.gamepad_radio, 'Value', ExpCfgParams.gamepad_radio);
set(handles.keyboard_radio, 'Value', ExpCfgParams.keyboard_radio);

set(handles.stimpath, 'String',[pwd,'\tempStimulus\']);
set(handles.ok_button, 'Enable', 'on');

if get(handles.quest_radio, 'Value') == 1;
    set(handles.priorSD, 'Visible', 'on');
    set(handles.text14, 'Visible', 'on', 'String', 'Beta');
    set(handles.beta, 'Visible', 'on');
    set(handles.text16, 'Visible', 'on');
    set(handles.delta, 'Visible', 'on');
    set(handles.text19, 'Visible', 'on');
    set(handles.pCorrect, 'Visible', 'on');
    set(handles.text17, 'String', 'Prior SD', 'Visible', 'on');
    set(handles.text13, 'Visible', 'on');
    set(handles.npresent, 'Visible', 'on', 'Enable', 'on');
    set(handles.quest_params, 'Title', 'QUEST Parameters');
    set(handles.text14, 'Visible', 'on', 'String', 'Beta');
    set(handles.beta, 'Visible', 'on', 'String', 3.5);
    set(handles.text16, 'Visible', 'on', 'String', 'Delta');
%     set(handles.delta, 'Visible', 'on', 'String', 0.01);
end

if get(handles.interleave_check, 'Value')==1;
    set(handles.low_pCorrect, 'Visible', 'on');
    set(handles.high_pCorrect, 'Visible', 'on');
else
    set(handles.low_pCorrect, 'Visible', 'off');
    set(handles.high_pCorrect, 'Visible', 'off');
end

if get(handles.gamepad_radio, 'Value')==0;
    %change text color
    set(handles.next_text, 'ForegroundColor', [1 1 1]);
    set(handles.yes_text, 'ForegroundColor', [1 1 1]);
    set(handles.repeat_text, 'ForegroundColor', [1 1 1]);
    set(handles.no_text, 'ForegroundColor', [1 1 1]);
    %enable text boxes
    set(handles.kb_next, 'Enable', 'on');
    set(handles.kb_yes, 'Enable', 'on');
    set(handles.kb_repeat, 'Enable', 'on');
    set(handles.kb_no, 'Enable', 'on');
else
    %change text color
    set(handles.next_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.yes_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.repeat_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.no_text, 'ForegroundColor', [0.5 0.5 0.5]);
    %disable text boxes
    set(handles.kb_next, 'Enable', 'off');
    set(handles.kb_yes, 'Enable', 'off');
    set(handles.kb_repeat, 'Enable', 'off');
    set(handles.kb_no, 'Enable', 'off');
end

%make 2 afc choice inactive for now;
set(handles.two_afc_radio, 'Enable', 'inactive');

set(handles.vidprefix, 'String', VideoParams.vidprefix);
user_entry = get(handles.initials,'String');
if get(handles.auto_prefix, 'Value') == 1
    set(handles.vidprefix, 'String', user_entry)
elseif get(handles.auto_prefix, 'Value') == 0
    %do nothing
end


% --- Executes when selected object is changed in stim_color_panel.
function stim_color_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in stim_color_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.red_stim_radio, 'Value')==1;
    set(handles.green_stim_radio, 'Value', 0);
elseif get(handles.green_stim_radio, 'Value') ==1;
    set(handles.red_stim_radio, 'Value', 0);
else
end


% function numrepeat_Callback(hObject, eventdata, handles)
% hObject    handle to numrepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numrepeat as text
%        str2double(get(hObject,'String')) returns contents of numrepeat as a double


% % --- Executes during object creation, after setting all properties.
% function numrepeat_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to numrepeat (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes during object creation, after setting all properties.
function stim_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in stim_panel.
function stim_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in stim_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)



function kb_next_Callback(hObject, eventdata, handles)
% hObject    handle to kb_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kb_next as text
%        str2double(get(hObject,'String')) returns contents of kb_next as a double


% --- Executes during object creation, after setting all properties.
function kb_next_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kb_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kb_yes_Callback(hObject, eventdata, handles)
% hObject    handle to kb_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kb_yes as text
%        str2double(get(hObject,'String')) returns contents of kb_yes as a double


% --- Executes during object creation, after setting all properties.
function kb_yes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kb_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kb_no_Callback(hObject, eventdata, handles)
% hObject    handle to kb_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kb_no as text
%        str2double(get(hObject,'String')) returns contents of kb_no as a double


% --- Executes during object creation, after setting all properties.
function kb_no_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kb_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function kb_repeat_Callback(hObject, eventdata, handles)
% hObject    handle to kb_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kb_repeat as text
%        str2double(get(hObject,'String')) returns contents of kb_repeat as a double


% --- Executes during object creation, after setting all properties.
function kb_repeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kb_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupilsize_Callback(hObject, eventdata, handles)
% hObject    handle to pupilsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupilsize as text
%        str2double(get(hObject,'String')) returns contents of pupilsize as a double



function fieldsize_Callback(hObject, eventdata, handles)
% hObject    handle to fieldsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldsize as text
%        str2double(get(hObject,'String')) returns contents of fieldsize as a double



function videodur_Callback(hObject, eventdata, handles)
% hObject    handle to videodur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of videodur as text
%        str2double(get(hObject,'String')) returns contents of videodur as a double



function presentdur_Callback(hObject, eventdata, handles)
% hObject    handle to presentdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of presentdur as text
%        str2double(get(hObject,'String')) returns contents of presentdur as a double
set(handles.presentdur_frames, 'String', num2str(round(str2double(get(handles.frame_rate, 'String')).*(str2double(get(handles.presentdur, 'String')))./1000)));


function iti_Callback(hObject, eventdata, handles)
% hObject    handle to iti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iti as text
%        str2double(get(hObject,'String')) returns contents of iti as a double


function vidprefix_Callback(hObject, eventdata, handles)
% hObject    handle to vidprefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vidprefix as text
%        str2double(get(hObject,'String')) returns contents of vidprefix as a double



function priorSD_Callback(hObject, eventdata, handles)
% hObject    handle to priorSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of priorSD as text
%        str2double(get(hObject,'String')) returns contents of priorSD as a double



function pCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pCorrect as text
%        str2double(get(hObject,'String')) returns contents of pCorrect as a double



function thresholdGuess_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdGuess as text
%        str2double(get(hObject,'String')) returns contents of thresholdGuess as a double



function beta_Callback(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta as text
%        str2double(get(hObject,'String')) returns contents of beta as a double



function delta_Callback(hObject, eventdata, handles)
% hObject    handle to delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta as text
%        str2double(get(hObject,'String')) returns contents of delta as a double



function npresent_Callback(hObject, eventdata, handles)
% hObject    handle to npresent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npresent as text
%        str2double(get(hObject,'String')) returns contents of npresent as a double



function gain_Callback(hObject, eventdata, handles)
% hObject    handle to gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gain as text
%        str2double(get(hObject,'String')) returns contents of gain as a double



function green_x_offset_Callback(hObject, eventdata, handles)
% hObject    handle to green_x_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_x_offset as text
%        str2double(get(hObject,'String')) returns contents of green_x_offset as a double



function green_y_offset_Callback(hObject, eventdata, handles)
% hObject    handle to green_y_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_y_offset as text
%        str2double(get(hObject,'String')) returns contents of green_y_offset as a double



function red_x_offset_Callback(hObject, eventdata, handles)
% hObject    handle to red_x_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_x_offset as text
%        str2double(get(hObject,'String')) returns contents of red_x_offset as a double



function red_y_offset_Callback(hObject, eventdata, handles)
% hObject    handle to red_y_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_y_offset as text
%        str2double(get(hObject,'String')) returns contents of red_y_offset as a double

function stim_midpoint_Callback(hObject, eventdata, handles)
% hObject    handle to stim_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_midpoint as text
%        str2double(get(hObject,'String')) returns contents of stim_midpoint as a double


% --- Executes during object creation, after setting all properties.
function stim_midpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_stims_riccos_Callback(hObject, eventdata, handles)
% hObject    handle to num_stims_riccos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_stims_riccos as text
%        str2double(get(hObject,'String')) returns contents of num_stims_riccos as a double


% --- Executes during object creation, after setting all properties.
function num_stims_riccos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_stims_riccos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function presentdur_frames_Callback(hObject, eventdata, handles)
% hObject    handle to presentdur_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of presentdur_frames as text
%        str2double(get(hObject,'String')) returns contents of presentdur_frames as a double
set(handles.presentdur, 'String', num2str(1000*str2double(get(handles.presentdur_frames, 'String'))./str2double(get(handles.frame_rate, 'String'))));


% --- Executes during object creation, after setting all properties.
function presentdur_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to presentdur_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_rate_Callback(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_rate as text
%        str2double(get(hObject,'String')) returns contents of frame_rate as a double
set(handles.presentdur_frames, 'String', num2str(round(str2double(get(handles.frame_rate, 'String')).*(str2double(get(handles.presentdur, 'String')))./1000)));


% --- Executes during object creation, after setting all properties.
function frame_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interleave_check.
function interleave_check_Callback(hObject, eventdata, handles)
% hObject    handle to interleave_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interleave_check
if get(handles.interleave_check, 'Value')==1;
    set(handles.low_pCorrect, 'Visible', 'on');
    set(handles.high_pCorrect, 'Visible', 'on');
else
    set(handles.low_pCorrect, 'Visible', 'off');
    set(handles.high_pCorrect, 'Visible', 'off');
end


function high_pCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to high_pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_pCorrect as text
%        str2double(get(hObject,'String')) returns contents of high_pCorrect as a double


% --- Executes during object creation, after setting all properties.
function high_pCorrect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high_pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function low_pCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to low_pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_pCorrect as text
%        str2double(get(hObject,'String')) returns contents of low_pCorrect as a double


% --- Executes during object creation, after setting all properties.
function low_pCorrect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_pCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in device_panel.
function device_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in device_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.gamepad_radio, 'Value')==0;
    %change text color
    set(handles.next_text, 'ForegroundColor', [1 1 1]);
    set(handles.yes_text, 'ForegroundColor', [1 1 1]);
    set(handles.repeat_text, 'ForegroundColor', [1 1 1]);
    set(handles.no_text, 'ForegroundColor', [1 1 1]);
    %enable text boxes
    set(handles.kb_next, 'Enable', 'on');
    set(handles.kb_yes, 'Enable', 'on');
    set(handles.kb_repeat, 'Enable', 'on');
    set(handles.kb_no, 'Enable', 'on');
else
    %change text color
    set(handles.next_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.yes_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.repeat_text, 'ForegroundColor', [0.5 0.5 0.5]);
    set(handles.no_text, 'ForegroundColor', [0.5 0.5 0.5]);
    %disable text boxes
    set(handles.kb_next, 'Enable', 'off');
    set(handles.kb_yes, 'Enable', 'off');
    set(handles.kb_repeat, 'Enable', 'off');
    set(handles.kb_no, 'Enable', 'off');
end
    


% --- Executes when main_CFG_figure is resized.
function main_CFG_figure_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to main_CFG_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
