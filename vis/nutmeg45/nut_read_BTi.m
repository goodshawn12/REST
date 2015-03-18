function [epochs,f0,no_points,start_time]=nut_read_BTi(filename,trials);
% [EPOCHS,F0,NO_POINTS,START_TIME] = 
%	NUT_READ_BTI(FILENAME, [TRIALS])
%
% reads BTi-format MEG data
% input:
%   FILENAME -- duh
%   TRIALS -- optional vector containing desired trial numbers;
%           if unspecified, returns all trials
% output:
%   EPOCHS -- MEG data matrix (no_channels by no_points by no_trials)
%   F0 -- sampling frequency
%   START_TIME -- time at first data point
%
% First created: 06/04/03 by SSN and SSD

tic

fid = fopen(filename,'r','b');
no_channels=fread(fid,1,'float');
no_epochs=fread(fid,1,'float');
sampling_period=fread(fid,1,'float');
start_time=fread(fid,1,'float');

epoch_no=fread(fid,1,'float');
no_points=fread(fid,1,'float');

if(nargin<2)
    trials = 1:no_epochs;
end

no_trials = length(trials);

epochs = zeros(no_channels,no_points,no_trials);


% epochs = zeros(no_channels,no_points,no_epochs);

% epochs(:,:,1)=fread(fid,[no_channels,no_points],'float');
% for trialnum=2:(no_epochs)
% 	fseek(fid,8,0);
% 	epochs(:,:,trialnum)=fread(fid,[no_channels,no_points],'float');
% end

h = waitbar(0,'Please wait, loading big-ass file...');

for i=1:no_trials;
    skipbytes=(trials(i)-1)*(no_channels*no_points+2)*4+24;
	fseek(fid,skipbytes,-1);
    epochs(:,:,i)=fread(fid,[no_channels,no_points],'float');
    waitbar(i/no_trials,h);
end

close(h);

f0 = 1/sampling_period;
toc
