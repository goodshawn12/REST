function pop = fcm_align2freqpeak(pop,api,wid)
% FCM_ALIGN2FREQPEAK   calculates population data relative to a frequency bin.
%  pop = fcm_align2freqpeak(pop,fbi,width)
%
% FBI    index of frequency bin to which data should be aligned for each subject.
% WIDTH  band width in Hz before and after the selected frequency bin.
%        Default is 4 Hz

if ~isvector(api), help fcm_align2freqpeak, error('FBI must be a vector.'), end
if nargin<3, wid=4; end

fsr = pop.timepts(2)-pop.timepts(1);
[numsubj,numvox,numfrq]=size(pop.s);

if numsubj~=length(api), error('FBI must be of same length as number of subjects.'), end

nas=zeros(numsubj,numvox,2*(wid/fsr)+1);
for k=1:numsubj
    nas(k,:,:,:)=pop.s(k,:,api(k)-(wid/fsr):api(k)+(wid/fsr));
end

pop.s=nas;
pop.timepts = [-wid:fsr:wid]';
pop.timewindow = [pop.timepts-fsr pop.timepts+fsr];
