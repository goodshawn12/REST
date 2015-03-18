function nut_import_Cartool_EPH(pa)
% Imports EPH data and converts to NUTMEG format
%
%   nut_import_Cartool_EPH(folder)
%
% folder   name of folder containing all EPH files of a given subject and condition.

global nuts

fieldtripcode=which('ft_senslabel');
if isempty(fieldtripcode)
    error('Please download, and place in your Matlab path, Fieldtrip code');
end

[labels,baseline,system]=getinput;

if ~isempty(labels)
    goodlab=ft_channelselection('gui',labels); pause(.1)
    good = find(ismember(labels,goodlab)); clear goodlab
    %bad = setdiff(1:length(labels),good);
end

% List all files
D=dir(pa); D(1:2)=[];
D={D.name};

numtrials=length(D);

% Read Header of first file
hdr=textread(fullfile(pa,D{1}),'%d',3);

% Read Data of all files in folder
fprintf('Reading Data: ')
data = zeros(hdr(2),hdr(1),numtrials);
for k=1:numtrials
    data(:,:,k) = readeph(fullfile(pa,D{k}));  
    if ~rem(k,5), fprintf('.'), end
end
fprintf('\n')

% Remove bad channels
if ~isempty(labels)
    data=data(:,good,:);
else
    good=1:hdr(1);
    labels=cell(hdr(1),1);
    for k=1:hdr(1)
        labels{k}=sprintf('E%03d',k);
    end
end

% Output
nuts.meg = struct('data', data, 'latency', ([0:hdr(2)-1]' ./ hdr(3)).*1000 - baseline, ...
    'srate', hdr(3), 'goodchannels', good, ...
    'sensor_labels', {labels}, ...
    'filename', [strtok(D{1},'.') '.Epoch all.eph'], 'eegflag', true, ...
    'system', system);

aw=questdlg('Do you want to average reference your data?','NUTMEG question','Yes','No','Yes');
if strcmp(aw,'Yes')
    nuts.meg=nut_eegref(nuts.meg,'AVG');
end

%---------------
function out = readeph(fn)

fid = fopen(fn,'rt');
C = textscan(fid,'','headerlines',1);
fclose(fid);

out=cat(2,C{:});

%----------------
function [labels,baseline,system]=getinput

ftchoice={
    'biosemi64'
    'biosemi128'
    'biosemi256'
    'bti148'
    'bti148_planar'
    'bti248'
    'bti248_planar'
    'btiref'
    'ctf151'
    'ctf151_planar'
    'ctf275'
    'ctf275_planar'
    'ctfheadloc'
    'ctfref'
    'eeg1005'
    'eeg1010'
    'eeg1020'
    'egi128'
    'egi256'
    'egi32'
    'egi64'
    'ext1020'
    'neuromag122'
    'neuromag122alt'
    'neuromag306'
    'neuromag306alt'
    'itab153'
    'itab153_planar'
    'yokogawa160'
    'yokogawa160_planar'
    'electrode'};

global temp
temp=[];
temp.infig = figure('units','normalized','position',[.45 .35 .15 .13], ...
    'toolbar','none','menubar','none','numbertitle','off');
uicontrol(temp.infig,'Style','text','string','Electrode system', ...
    'units','normalized','position',[.05 .8 .35 .13],'horizontalalignment','left','backgroundcolor',[.8 .8 .8]);
uicontrol(temp.infig,'Style','popupmenu','string',ftchoice,'value',2,'tag','syspop', ...
    'units','normalized','position',[.05 .67 .9 .13],'horizontalalignment','left','backgroundcolor','w');
uicontrol(temp.infig,'Style','text','string','Time origin (0 ms)', ...
    'units','normalized','position',[.05 .4 .35 .13],'horizontalalignment','left','backgroundcolor',[.8 .8 .8]);
uicontrol(temp.infig,'Style','edit','units','normalized','tag','bsedit', ...
    'position',[.45 .4 .2 .13],'backgroundcolor','w');
uicontrol(temp.infig,'Style','text','string','ms after start', ...
    'units','normalized','position',[.7 .4 .3 .13],'horizontalalignment','left','backgroundcolor',[.8 .8 .8]);
uicontrol(temp.infig,'Style','pushbutton','string','Ok', ...
    'callback','global temp; handles=guihandles(temp.infig); temp.baseline=get(handles.bsedit,''string''); temp.system=get(handles.syspop,''value''); delete(gcf);', ...
    'units','normalized','position',[.3 .1 .4 .2]);
uiwait(temp.infig);
pause(.1)

system=ftchoice{temp.system};
labels=ft_senslabel(system);
baseline = str2num(temp.baseline);
if isempty(baseline), baseline=0; end
clear global temp




