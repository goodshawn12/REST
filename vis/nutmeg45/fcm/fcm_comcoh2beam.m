function fcm_comcoh2beam(comcohfile,wfile)
% FCM_COMCOH2BEAM  checks your current FCM configuration and creates output 
%  beam file from functional connectivity values.  
%
%    fcm_comcoh2beam(comcohfile,¦wfile¦)
%
%  COMCOHFILE  name of file containing functional connectivity results.
%  WFILE       (optional) File containing the spatial filter matrix W. If
%              this input is given, the function will try to create an
%              additional set of files that are corrected for artifacts due
%              to spatial leakage of the inverse solution. The correction
%              subtracts the spatial filter correlation matrix from the
%              data functional connectivity matrix.


global fuse nuts

if (~isstruct(fuse) || ~isfield(fuse,'seed'))
    error('You must set configuration first (with fcm_set_config or fcm_config_gui).')
end
if (~isstruct(nuts) || ~isfield(nuts,'voxels'))
    error('You must open the corresponding session file first.')
end
if nargin<2, wfile=[]; end

switch fuse.connection
    case 'All'
        switch fuse.seed
            case 'All'
                meanfun = 'fcm_comcohA2meanbeam';
                Lfun    = 'fcm_comcohA2Lbeam';
            case {'Selected' 'Selected+Contralateral'}
                meanfun = 'fcm_comcohS2meanbeam';
                if ~isempty(strmatch('L',fuse.output,'exact')) && ~strcmpi(fuse.seed,'Selected+Contralateral')
                    warning('fcm:noLimage','L-images can only be computed with seed voxels set to ''All'' or ''Selected+Contralateral''');
                    Lfun= [];
                else
                    Lfun= 'fcm_comcohSC2Lbeam';
                end            
            case 'Extracerebral'
                meanfun = 'fcm_comcohE2meanbeam';
                Lfun    = [];
                if ~isempty(strmatch('L',fuse.output,'exact'))
                    warning('fcm:noLimage','L-images can only be computed with seed voxels set to ''All'' or ''Selected+Contralateral''');
                end     
            otherwise
              error('Invalid seed voxel setting.')  
        end
        
    case 'Grid'
        meanfun = 'fcm_comcohG2meanbeam';
        if ~isempty(strmatch('L',fuse.output,'exact'))
            switch fuse.seed
                case 'Selected+Contralateral'
                    Lfun = 'fcm_comcohSC2Lbeam';
                case 'All'
                    Lfun = 'fcm_comcohAG2Lbeam';
                otherwise
                    Lfun = [];
                    warning('fcm:noLimage','L-images can only be computed with seed voxels set to ''All'' or ''Selected+Contralateral''');
            end
        else
            Lfun = [];
        end
                
    otherwise
        error('Invalid connection voxel setting.')
end

correctnout(meanfun,Lfun,comcohfile,wfile)

%-------------------
function correctnout(meanfun,Lfun,comcohfile,wfile)
global fuse

produceoutput(meanfun,Lfun,comcohfile)

% Correct spatial leakage
if ~isempty(wfile) && ~strcmpi(fuse.seed,'Extracerebral') && any(strcmp(fuse.funconn,{'ccohere' 'nccohere' 'ampcorr'}))   % Correction for PLI not necessary
    load(comcohfile);
    load(wfile);
    if ~exist('W','var')
        if exist('Wact','var'), W=Wact; 
        else warning('Not a valid weights file, cannot correct for spatial leakage.'), return
        end
    end
    if ndims(W)>2 && size(W,2)>1
        warning('No scalar filter weights provided, cannot correct for spatial leakage.'), return
    end    
    CC=fcm_correctspatialleakage(CC,W); clear W
    [pa,fi,dum]=fileparts(comcohfile);
    comcohfile = fullfile(pa,[fi 'C']);
    save(comcohfile,'CC'), clear CC
    produceoutput(meanfun,Lfun,comcohfile)
end
    
%-----------------------
function produceoutput(meanfun,Lfun,comcohfile)
global fuse nuts

[pa,fi,dum]=fileparts(comcohfile);
beamfile = fullfile(pa,['s_beamtf_' fi]);
if strcmpi(fuse.connectionavg,'Seed')
    beamfile=[beamfile '_seed'];
    if isfield(nuts,'selvox') && isfield(nuts.selvox,'ipsilab')
        if iscell(nuts.selvox.ipsilab), lab=[nuts.selvox.ipsilab{:}];
        else lab=nuts.selvox.ipsilab; 
        end
        beamfile=[beamfile strrep(lab,'_','')];
    end
elseif ~strcmpi(fuse.connectionavg,'all')
    beamfile = [beamfile '_' lower(fuse.connectionavg)];
end
if ~isempty(strmatch('Mean',fuse.output,'exact'))
    if strcmpi(fuse.seed,'Extracerebral') || strcmpi(fuse.connection,'Grid')
        beam = feval(meanfun,comcohfile);
    else
        beam = feval(meanfun,comcohfile,fuse.connectionavg);
    end
    if length(beam)>1
        BEAM=beam; clear beam
        for k=1:length(BEAM);
            beam = BEAM(k);
            save([beamfile '_ext' int2str(k)],'beam')
            if ~isempty(strmatch('Z',fuse.output,'exact')) 
                beam = fcm_znormbeam(beam);
                save([beamfile '_ext' int2str(k) 'Z'],'beam')
            end            
        end
    else
        save(beamfile,'beam')	
        if isfield(beam,'baseline')
            bb=beam;
            beam=fcm_basebeam(beam);
            if length(beam)>1
                txt = {'imag' 'mag2' 'real'};
                BEAM=beam;
                for k=1:length(beam)
                    beam = BEAM(k);
                    save([beamfile '_' txt{k} 'base'],'beam')
                end
            else
                save([beamfile '_base'],'beam')
            end
            beam=bb;
        end        
        if ~isempty(strmatch('Z',fuse.output,'exact')) 
            beam = fcm_znormbeam(beam);
            save([beamfile 'Z'],'beam')
            if isfield(beam,'baseline')
                beam=fcm_basebeam(beam);
                if length(beam)>1
                    txt = {'imag' 'mag2' 'real'};
                    BEAM=beam;
                    for k=1:length(beam)
                        beam = BEAM(k);
                        save([beamfile 'Z_' txt{k} 'base'],'beam')
                    end
                else
                    save([beamfile 'Z_base'],'beam')
                end                
            end
        end
    end
end
if ~isempty(strmatch('L',fuse.output,'exact')) && ~isempty(Lfun)
    beam = feval(Lfun,comcohfile);
    save([beamfile 'L'],'beam')
end
