function nut_beam2cartool(beam,risfilename,freqband)
% NUT_BEAM2CARTOOL  exports NUTMEG beam structure to Cartool for visualization.
%
%     nut_beam2cartool(beam,risfilename,frequencyband)
%
% All input arguments are optional.
%
% beam          NUTMEG beam structure. If no input argument is given, the data 
%               currently displayed in nut_results_viewer is exported.
% risfilename   name of RIS file to be saved. Default is GUI input.
% frequencyband number of frequency band to be exported. Default is the
%               band currently selected in nut_results_viewer or 1.

global rivets

% Check inputs
if nargin<1 || isempty(beam)            % if beam input not given
    %clear beam                         % this is necessary to avoid MATLAB warning.
    global beam
    if ~isfield(rivets,'s')
        error('You must either have the results displayed in the NUTMEG viewer or specify a beam structure as first input argument.')
    end
end
if isfield(beam,'sa')        % legacy compatibility
    beam.s{1}=beam.sa{1};
    if( isfield(beam,'sc') && any(beam.sc{1}(:)~=1) )
        beam.s{2}=beam.sc{1};
    end
    if isfield(beam,'n')
        beam.s{3}=beam.n{1}; 
    end
end

if (nargin<2 || isempty(risfilename))
    [fRIS, pRIS] = uiputfile([rivets.beamfilename '.ris'], 'Export as...');
    if isequal(fRIS,0), return, end
    risfilename=fullfile(pRIS,fRIS);
end

if nargin<3
    if isfield(rivets,'freqselect')     % if nut_timef_viewer open and displaying data
        freqband = rivets.freqselect;
    else
        freqband=1;
        if size(beam.s{1},3)>1
            warning('Exporting first frequency band only!')
        end
    end
end

% Prepare Data
if isfield(rivets,'s')      % if nut_timef_viewer open and displaying data
   data = rivets.s(:,:,freqband); 
else    
    if length(beam.s)>1
        clear global stats
        beam = nut_calc_beamratio(beam);
    end
    data = beam.s{1}(:,:,freqband);
end

if isfield(rivets,'threshold')  % Statistical thresholding
    notsignificant = ~any(rivets.threshold,3);
    data(notsignificant) = 0;
end

% Write to file
nut_write2cartoolris(risfilename,data,beam.srate);

% Export solution points if requested
answer=questdlg('Would you like to export the voxel coordinates?','NUTMEG question','Yes','No','No');
switch answer
    case 'Yes'
        answer2=questdlg('To which MRI orientation?','NUTMEG question','SPM','SMAC','origin','SPM');
        
        [fSPI, pSPI] = uiputfile([rivets.beamfilename '.spi'], 'Export as...');
        if isequal(fSPI,0), return, end
        spifilename=fullfile(pSPI,fSPI);
        nut_write2cartoolspi(spifilename,beam,answer2);
end
