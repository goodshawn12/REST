function beam=nut_oldbeam2newbeam(beam)
% converts nutmeg 2.0 files to nutmeg 3.0 files


if(isfield(beam,'tf'))  
    beam.s{1}=beam.tf;
    beam = rmfield(beam,'tf');
end
if isfield(beam,'s')
    % s_beamtf used to have this format
    if(all(size(beam.s)==[1 3]))
        beam.sa{1}=beam.s{1};
        beam.sc{1}=beam.s{2};
        beam.n{1}=beam.s{3};

        % sarang's evoked s_beam had this format (when there are 2 conditions)
    elseif(all(size(beam.s)==[1 2]))
        beam.sa=beam.s
        beam.timepts = beam.timewindow;
        beam = rmfield(beam,'timewindow');
        
    % This is the format of the CTF style SAM and of the fcm stuff   
    elseif(all(size(beam.s)==[1 1]))
        beam.sa{1}=beam.s{1};
        beam.sc{1}=ones(size(beam.sa{1}),'single');
        
    else
        error('your weird format not yet supported -- hey, why don''t you code it up? :)')
    end
    beam=rmfield(beam,'s');
elseif isfield(beam,'s_th') % compatible with nut_activation_viewer
    beam.timepts = beam.timewindow;
    beam = rmfield(beam,'timewindow');
    beam.sa{1}(:,:,1,1)=beam.s_th;
    if ~all(beam.s_ph(:)==0)
        beam.sa{1}(:,:,1,2)=beam.s_ph;
    end
    beam = rmfield(beam,'s_th');
    beam = rmfield(beam,'s_ph');
    if isfield(beam,'s_z')
        if ~all(beam.s_z(:)==0)
            beam.sa{1}(:,:,1,3)=beam.s_z;
        end
        beam = rmfield(beam,'s_z');
    end
else
    error('your weird format not yet supported -- hey, why don''t you code it up? :)')
end

if(isfield(beam,'p'))
    beam.normcdf.p_uncorr = beam.p;
    beam.normcdf.tail = 'both';
    beam.normcdf.cutoff = 0.01;
    beam = rmfield(beam,'p');
end


% evoked results did not previously have bands field, now required
if(~isfield(beam,'bands'))
    beam.bands = [0 inf];
end

% This seems to be required now as well...
if ~isfield(beam,'params') || ~isfield(beam.params,'beamformertype')
    beam.params.beamformertype='unknown';
end

% stats legacy compatibility
isask=false;
if isfield(beam,'snpm') && ~isequal(size(beam.snpm.T,1),size(beam.sa{1},1))
    beam.snpm.T=reshape(beam.snpm.T',[size(beam.s{1},1) 1 size(beam.sa{1},3)]);
    if isfield(beam.snpm,'p_uncorr_pos'), beam.snpm.p_uncorr_pos=reshape(beam.snpm.p_uncorr_pos',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]); end
    if isfield(beam.snpm,'p_uncorr_neg'), beam.snpm.p_uncorr_neg=reshape(beam.snpm.p_uncorr_neg',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]); end
    if isfield(beam.snpm,'p_corr_pos'),   beam.snpm.p_corr_pos=reshape(beam.snpm.p_corr_pos',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]); end
    if isfield(beam.snpm,'p_corr_neg'),   beam.snpm.p_corr_neg=reshape(beam.snpm.p_corr_neg',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]); end
    isask=true;
end
if isfield(beam,'ttestFDR') && ~isequal(size(beam.ttestFDR.T,1),size(beam.sa{1},1))
    beam.ttestFDR.T=reshape(beam.ttestFDR.T',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]);
    beam.ttestFDR.p_uncorr=reshape(beam.ttestFDR.p_uncorr',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]);
    if isfield(beam.ttestFDR.p_corr), beam.ttestFDR.p_corr=reshape(beam.ttestFDR.p_corr',[size(beam.sa{1},1) 1 size(beam.sa{1},3)]); end
    isask=true;
end
if isask
    warndlg('Your statistics are in an old format, which will not be supported much longer...!') 
end

