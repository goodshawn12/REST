function nuts = nut_bst2nuts(nuts,forwfile,chanfile,voxelstruct)
% converts BrainStorm's data format into nutmeg session structure
% function nuts = nut_bst2nuts(megfile, forwfile,covfile)
% megfile = file containing MEG data (FIFF format)
% forwfile = file containing forwards solution, i.e., the lead field
% covfile = file containing covariance (but more importantly, also
% containing source locations for forwfile)
% calls Brainstorm's read_gain function, provided below under GNU GPL


% global nuts
oldbrainstorm=0;

tmp = load(forwfile);
forw = tmp.forw_LCMVmethod;

load_data = 0;
if(load_data)
    [pathname,filebase,fileext]=fileparts(megfile);
    eventfile = fullfile(pathname,[filebase '-eve' fileext]);

    megdata_filename = megfile;
    nuts.meg.filename = megfile;

    mnecode=which('fiff_setup_read_raw');
    if isempty(mnecode)
        error('Please place the MNE Matlab software in your Matlab path. If you need to download it, please see http://nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php');
    end

    [ fid, tree ] = fiff_open(megfile);
    [ info, meas ] = fiff_read_meas_info(fid,tree);
    FIFF = fiff_define_constants();
    raw = fiff_dir_tree_find(meas,FIFF.FIFFB_RAW_DATA);

    if isempty(raw)  % for averaged data *.fif file
        avdata=fiff_read_evoked_all(megdata_filename);
        if length(avdata.evoked)>1
            whichepoch = input(['Which epoch to use of the ' num2str(length(avdata.evoked)) ' available?  ']);
        end
        megtmp=avdata.evoked(whichepoch).epochs';
        latency=avdata.evoked(whichepoch).times;
        [info]=fiff_read_meas_info(megdata_filename);
        %             sensor_labels=info.ch_names;
        jj=1;
        for ii=1:length(info.chs)
            if info.chs(ii).kind==1 % MEG channel type
                coil_coord(jj,:)=info.chs(ii).loc;
                meg(:,jj)=megtmp(:,ii);
                sensor_labels{jj}=info.chs(ii).ch_name;
                jj=jj+1;
            end
        end

    else % for raw *.fif data
        blah=fiff_setup_read_raw(megdata_filename);
        [meg,latency] = fiff_read_raw_segment(blah);

        jj=1;
        for ii=1:length(info.chs)
            if info.chs(ii).kind==1 % MEG channel type
                coil_coord(jj,:)=info.chs(ii).loc;
                ind(jj)=ii;
                %                     meg(:,jj)=megtmp(:,ii);
                sensor_labels{jj}=info.chs(ii).ch_name;
                jj=jj+1;
            end
        end
        kill = setdiff(1:length(info.chs),ind);
        meg(kill,:)=[];
        meg = meg';
    end
    normal_direction=zeros(size(coil_coord,1),1);  % for now this will be meaningless
    coil_coord2=zeros(size(coil_coord));
    srate=info.sfreq;
    nuts.meg.sensorCoord = coil_coord;
    nuts.meg.sensor_labels = sensor_labels;
    nuts.meg.srate = srate;
    % nuts.meg.sensorCoord = data.grad.pnt*1000;
    nuts.meg.Gcoef = 0;
    % BTi data in meters, but other data formats might have different units???
    % this could be serious problem with fieldtrip's leadfield computation
    % probably *does not* scale!

    % data.grad.tra contains synthetic third gradiometer coefficients
    nuts.meg.grad_order = 0;  % for BTi data
    % nuts.meg.refSensorCoord;
    % nuts.meg.refSensorOrient;
    % nuts.meg.ref_sensor_labels;

    % nuts.meg.sensorOrient = data.grad.ori;

    nuts.meg.latency = latency';
    nuts.meg.data = meg;
    clear meg
    nuts.meg.goodchannels = [1:size(nuts.meg.data,2)];

    nuts.coreg.mripath = which('blank.img');
    nuts.coreg.meg2mri_tfm = eye(4);

    markerdelay = 0.038 * nuts.meg.srate   %%%% NOTE: this delay is only for NeuroSpin's projection system as of May 28, 2009!!!!
    eventlist = mne_read_events(eventfile);
    markerlist = unique(eventlist(:,3));
    for ii=1:length(markerlist)
        markerlatencies{ii}=eventlist(find(eventlist(:,3)==markerlist(ii)),1) + markerdelay;
    end
    nuts.meg.markersnmag.latencies = markerlatencies;
    nuts.meg.markersnmag.codes = markerlist;

end

% nuts.voxels = nut_nmag2meg(1000*Vertices{1}');
%nuts.voxels = nut_mri2meg(1000*Vertices{1}');


% FIXME: use nuts.meg.goodchannels here instead of 'kill'
% load(forwfile);
% calling Brainstorm's read_gain function, provided below under GNU GPL
if(oldbrainstorm)
    load(voxelsfile);
% FIXME: figure out voxelsize for real?
    nuts.voxelsize = [5 5 5];
    nuts.voxels = nut_mri2meg(1000*Vertices{1}');
    g = read_gain(forwfile)
    % g(kill,:)=[];
    switch(size(g,2))
        case size(nuts.voxels,1)
            nuts.Lp(:,1,:) = g;
        case 3*size(nuts.voxels,1)
            nuts.Lp(:,1,:) = g(:,2:3:end);
            nuts.Lp(:,2,:) = g(:,1:3:end);
            nuts.Lp(:,3,:) = g(:,3:3:end);
        otherwise
            error('I don''t get this crap.');
    end
else
    chans=load(chanfile);
    magnetometer_idx=strmatch('MEG MAG',{chans.Channel.Type});
    nuts.voxels = 1000*forw.GridLoc';
    nuts.voxelsize = [5 5 5] % it's not equally spaced but whatever
    switch(size(forw.Gain,2))
        case size(nuts.voxels,1)
            nuts.Lp(:,1,:) = forw.Gain;
        case 3*size(nuts.voxels,1)
            nuts.Lp(:,1,:) = forw.Gain(:,2:3:end);
            nuts.Lp(:,2,:) = forw.Gain(:,1:3:end);
            nuts.Lp(:,3,:) = forw.Gain(:,3:3:end);
        otherwise
            error('I don''t get this crap.');
    end
    
    nuts.Lp=nuts.Lp(magnetometer_idx,:,:);
end

    

return





function G = read_gain(FILENAME,chanID,srcID)
% READ_GAIN: Extract parts or a compete gain matrix from .bin binary gain file.
%
% USAGE:  G = read_gain(FILENAME,chanID,srcID);
%
% INPUT:
%    - FILENAME: A character string describing the name of the file containing the matrix 
%                The gain matrix file must be in the BrainStorm binary file format.
%    - chanID  : Optional vector of indices linked to the ROWS of the matrix (ie channels) to be extracted 
%                If chanID is LEFT EMPTY, the matrix is extracted for all channels
%    - srcID   : Optional vector of indices linked to the COLUMNS of the matrix to be extracted 
%                If srcID is LEFT EMPTY, matrix is extracted for all sources
% OUTPUT:
%    - G       : an array containing the forward fields of the gain (sub)matrix.
%
% SEE ALSO: GET_GAIN, LOAD_RAW, SAVE_RAW

% @=============================================================================
% This software is part of The BrainStorm Toolbox
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2008 BrainStorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm licence" at command prompt.
% =============================================================================@
%
% Authors: Sylvain Baillet, October 2002
% ----------------------------- Script History ---------------------------------
% SB     Oct-2000  Creation
% JCM 27-May-2004  Comments updating
% ------------------------------------------------------------------------------


% Open the file
fid = fopen(FILENAME ,'r','ieee-be'); 
% If file not found: add STUDIES directory
if fid < 0
    ProtocolInfo = bst_getContext('ProtocolInfo');
    fid = fopen(fullfile(ProtocolInfo.STUDIES, FILENAME) ,'r','ieee-be'); 
end

if fid < 0 
   error('HeadModel file not found.');
end

if nargin == 1
   G = get_gain(fid);
elseif nargin == 2
   G = get_gain(fid,srcID);
elseif nargin == 3
   if isempty(srcID) % ChanID is defined but srcID is left blank
      G = get_gain(fid);
   else
      G = get_gain(fid,srcID);
   end
   if ~isempty(chanID)
       G = G(chanID,:); % Keep only channels of interest
   end
end

fclose(fid);


  

