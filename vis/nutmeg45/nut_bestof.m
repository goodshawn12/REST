function nut_bestof(maxtrial,lofrq,hifrq)
% NUT_BESTOF  keeps trials showing most activity within a frequency range.
%
%   nut_bestof(num_trials_to_keep,lowfreq,highfreq)
%     all input arguments are optional

global nuts

if nargin<1 || isempty(maxtrial), maxtrial=300; end
if nargin<2 || isempty(lofrq), lofrq=6; end
if nargin<3, hifrq=14; end

if size(nuts.meg.data,3)<=maxtrial, return, end

idx = frqscore(nuts.meg.data,[lofrq hifrq]);
nuts.meg.data = nuts.meg.data(:,:,sort(idx(1:maxtrial)));

function idx=frqscore(data,frqlim)
% SCORE  calculates an index of frequency activity of data segments.
%    idx = frqscore(nuts.meg.data)

F=fft(data,256,1);
f=0:2:254;
F(129:end,:,:)=[];
F=abs(F).^2;
v=find(f>=frqlim(1) & f<=frqlim(2));
%b=find(f>19 & f<40);
b=setdiff(find(f<60),v);
score=mean(mean(F(v,:,:),1),2) ./ mean(mean(F(b,:,:),1),2);
score=squeeze(score);
[score,idx]=sort(score,1,'descend');
