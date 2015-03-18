function varargout = nut_about(varargin)
% function varargout = NUT_ABOUT(varargin)
% The function displays the 'ABOUT' dialog, which is
% described in the NUT_ABOUT.FIG file.

warning('off')  % hack so warning for "created by matlab 7.0" is turned off
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_about_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_about_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
warning('on')  % we want to see all other warnings.


%%%-------------------------
function nut_about_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_about is made visible.
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Color',[1 1 1]);

axes(handles.nutmeglogo_axes);
banner = imread(['nutmeg.png']); % Read the image file
info = imfinfo(['nutmeg.png']); % Determine the size of the image file
image(banner);
set(handles.nutmeglogo_axes,'Visible', 'off');

axes(handles.UCSFlogo_axes);
banner = imread(['UCSFlogo.png']); % Read the image file
info = imfinfo(['UCSFlogo.png']); % Determine the size of the image file
image(banner);
set([handles.nutmeglogo_axes handles.UCSFlogo_axes],'Visible', 'off');

viewcredits_button_Callback(hObject,eventdata,handles);  % load up with credits displayed




%%%-------------------------
function varargout = nut_about_OutputFcn(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
varargout{1} = handles.output;



% --- Executes on button press in viewlicense_button.
function viewlicense_button_Callback(hObject, eventdata, handles)
% hObject    handle to viewlicense_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nutmegpath = fileparts(which('nutmeg'));
fid=fopen([nutmegpath '/LICENSE.txt']);
linenum = 1;
while 1
    license{linenum}=fgetl(fid);
    if ~ischar(license{linenum})
        break
    end
    linenum = linenum + 1;
end
fclose(fid);
license = {license{1:(end-1)}};
set(handles.license_text,'FontSize',10,'FontName','Helvetica');
set(handles.license_text,'String',textwrap(handles.license_text,license));

set(handles.license_text,'Visible','On');
set(handles.viewlicense_button,'Visible','Off');
set(handles.viewcredits_button,'Visible','On');







% --- Executes on button press in viewcredits_button.
function viewcredits_button_Callback(hObject, eventdata, handles)
% hObject    handle to viewcredits_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if(strcmp(version('-release'),'14'))
    oumlat = char(246);
else
    oumlat = 'o';  % R13 doesn''t like the umlat, apparently
end

header = {'',...
'  Neurodynamic Utility Toolbox for Magnetoencephalography - NUTMEG',...
'  http://nutmeg.berkeley.edu',...
'',...
'  Developed by the UCSF Biomagnetic Imaging Laboratory and friends',...
'  http://bil.ucsf.edu/'};

authors = {'  Current Authors:',...
'     Sarang S. Dalal',...
'     Johanna M. Zumer',...
'     Adrian G. Guggisberg',...
'     Julia P. Owen',...
'     Kensuke Sekihara',...
'     Hagai T. Attias',...
'     Srikantan S. Nagarajan'};

formerauthors = {'','   Former Contributors:',...
'     Alex Wade',...
'     Vineet Agrawal',...
'     DoSik Hwang',...
'     Kenneth Hild'};

thanks = {'  Special Thanks:',...
'     Dave P. Wipf',...
'     Dave J. McGonigle',...
'     Anne M. Findlay',...
'     Darren Weber',...
'     John F. Houde',...
'     Tony Norcia',...
'     Janine Lupo',...
['     Gr' oumlat 'lsch'],...
'     Casa Noble',...
'     Yancy'};
footer = {['   ' char(169) '2005-2009 The Authors and the UCSF Biomagnetic Imaging Laboratory'],...
'   All Rights Reserved.'};
set([handles.authors_text handles.header_text handles.thanks_text],'FontSize',12,'FontName','Helvetica');
set(handles.authors_text,'String',textwrap(handles.authors_text,[authors formerauthors]));
set(handles.header_text,'String',textwrap(handles.header_text,header));
set(handles.thanks_text,'String',textwrap(handles.thanks_text,thanks));
set(handles.footer_text,'String',textwrap(handles.footer_text,footer));
set(handles.license_text,'Visible','Off');
set(handles.viewlicense_button,'Visible','On');
set(handles.viewcredits_button,'Visible','Off');

