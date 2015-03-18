function nut_write2ctf(ctffolder,data,chantype)
% NUT_WRITE2CTF: Replaces channels in CTF dataset with data in a MATLAB variable.
%
%    nut_write2ctf(ctffolder,meg_data,channel_type);
%
%    ctffolder      CTF dataset (*.ds)
%    data           Variable containing new data (time*channels*trials). The number 
%                   of samples, channels, and trials must match the numbers in the 
%                   original CTF dataset.
%    channel_type   Type of channels in 'meg_data' (optional, default is 'meg', 
%                   type 'help ctf_channel_select' for other options).

if nargin<3, chantype='meg'; end

[folderPath,folderName,folderExt] = fileparts(ctffolder);
% meg4file=fullfile(ctffolder,[folderName '.meg4'])

% read res4 file and find MEG channels
ctf  = ctf_read_res4(ctffolder);
if ~isequal(ctf.setup.number_samples,size(data,1))
    error('NUT_WRITE2CTF: data must have same number of samples as ctf dataset.')
end
if ~isequal(ctf.setup.number_trials,size(data,3))
    error('NUT_WRITE2CTF: data must have same number of trials as ctf dataset.')
end
CHAN = ctf_channel_select(ctf,chantype);
if ~isequal(length(CHAN),size(data,2))
    error(sprintf('NUT_WRITE2CTF: data must have same number of %s channels as ctf dataset.'),chantype)
end

% open meg4 file
cd(ctffolder);
[fid,msg]=fopen([folderName '.meg4'],'r+','s');
if fid < 0, error(msg); end


% replace MEG data in meg4 file
fprintf('\nNUT_WRITE2CTF: replacing %s channels in %s...\n',chantype,folderName)
header_bytes = 8;
trial_bytes = 4 * ctf.setup.number_samples * ctf.setup.number_channels;

for trial = 1:ctf.setup.number_trials

  % calculate the byte offset in the file for this trial
  trial_offset = header_bytes + ( trial_bytes * ( trial - 1 ) );
  
  for channel = 1:length(CHAN),
    channel_number = CHAN(channel);
    
    % calculate the channel offset in the current trial
    channel_offset = 4 * ctf.setup.number_samples * ( channel_number - 1 );
    
    % seek to the trial offset, relative to the beginning of the file
    fseek(fid,trial_offset,-1);
    
    % now seek to the channel offset, in the current trial
    fseek(fid,channel_offset,0);
    
    % write
    fwrite(fid,data(:,channel,trial),'int32');
  end

end
fclose(fid);
cd ..