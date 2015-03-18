function nut_importmat(data,srate,latency,lsc,sensorCoord,sensorOrient,sensor_labels)
% nut_importmat(matfile)
% nut_importmat(data,srate,lsc,sensorCoord)
% nut_importmat(data,srate,lsc,sensorCoord,sensor_labels)
%
% Imports MEG time series from a mat file or variable

global nuts

if(nargin==1)
    nuts.meg.filename = data;
    meg = load(data);
    if(all(ismember({'data','lsc','sensorCoord','sensorOrient','srate','latency'},fieldnames(meg))))
        nuts.meg = meg;
        disp('Good MEGer.')
    else
        prompt={'MEG Time Series, time x channels:'
                'Sampling Rate (Hz):'
                'Sample latencies (ms), time x 1:'
                'Local Sphere Center (mm):'
                'Sensor Locations (mm), channels x 3 x 2:'
                'Sensor Orientations (mm), channels x 3:'
                'Sensor Names (optional):'};
        name='Import MAT file';
        numlines=1;
        defaultanswer={ 'data'
                        'srate'
                        'latency'
                        'lsc'
                        'sensorCoord'
                        'sensorOrient'
                        'sensor_labels'   };
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        nuts.meg.data = getfield(meg,answer{1});
        nuts.meg.srate = getfield(meg,answer{2});
        nuts.meg.latency = getfield(meg,answer{3});
        nuts.meg.lsc = getfield(meg,answer{4});
        nuts.meg.sensorCoord = getfield(meg,answer{5});
        nuts.meg.sensorCoord = getfield(meg,answer{6});
        if(~isempty(answer{7}))
            nuts.meg.sensor_labels = getfield(meg,answer{8});
        end

        if(isempty(answer))
            return; % no input, let's blow this taco stand
        end
    end
else
    nuts.meg.filename = inputname(1);
    nuts.meg.data = data;
    nuts.meg.lsc = lsc;
    nuts.meg.sensorCoord = sensorCoord;
    nuts.meg.sensorOrient = sensorOrient;
    nuts.meg.latency = latency;
end

if(exist('sensor_labels','var'))
    nuts.meg.sensor_labels = sensor_labels;
else
    leftchans = find(nuts.meg.sensorCoord(:,2,1) >= 0);
    rightchans = find(nuts.meg.sensorCoord(:,2,1) < 0);
    for i = 1:length(leftchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string];
    end
    for i = 1:length(rightchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string];
    end
end

nuts.meg.goodchannels=1:size(nuts.meg.sensorCoord,1);

nut_enabler;
