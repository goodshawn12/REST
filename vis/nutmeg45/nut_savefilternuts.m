function nut_savefilternuts(sessionfile,lowfreq,hifreq,varargin)
% function nut_savefilternuts(sessionfile,lowfreq,hifreq,varargin)
%
% tfcov(sessionfile,paramfile,lowfreq,highfreq)
% wrapper script for nut_tfcov, intended for compilation into standalone
% varargin{1} is .mat file with filt. and params. stuff set in it

tic

if(ischar(lowfreq))
    lowfreq = str2num(lowfreq);
end
if(ischar(hifreq))
    hifreq = str2num(hifreq);
end

if(nargin<4)
    filt.class='butter';
    filt.order=4;
    filt.type='bp';
    filtname='butterbp4';
%     algo={'SAM'};
%     params.cn=1;
%     params.qsub=0;
else
    load(varargin{1}); % load additional params (like filter settings)
    [filtpath,filtname,ext] = fileparts(varargin{1});
%     algo=cell(length(varargin)-1,1);
%     [algo{1:size(algo,1)}]=varargin{2:end};
% algo={'SAM'};
end

nuts=load(sessionfile); % must contain nuts.Lp, nuts.voxels, nuts.voxelsize, nuts.coreg, nuts.meg
% backwards compatibility
if isfield(nuts,'nuts')
    nuts=nuts.nuts;
end
% load(paramfile);  % must contain time windows ("active" and "control")

% print run parameters (for preservation in log files)
disp(sessionfile)
% disp(paramfile)
filt

tic
for ii=1:size(nuts.meg.data,3)
    nuts.meg.data(:,:,ii) = nut_filter2(squeeze(nuts.meg.data(:,:,ii)),filt.class,filt.type,filt.order,lowfreq,hifreq,nuts.meg.srate,1);
end
% nuts.meg.data = nut_filter2(nuts.meg.data,filt.class,filt.type,filt.order,lowfreq,hifreq,nuts.meg.srate,1);
toc

[filepath,filename,ext] = fileparts(sessionfile);
outname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz.mat'];

meg=nuts.meg;
clear nuts
% nuts.meg.data=single(nuts.meg.data);
meg.data=single(meg.data);
stemp=whos('meg');
if stemp.bytes > 2e9
    save(outname,'-v7.3','meg');
else
    save(outname,'meg');
end

