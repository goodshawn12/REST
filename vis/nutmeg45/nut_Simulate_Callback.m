function nut_Simulate_Callback(answer)
% makes up some data with known sources to test beamforming routines
% borrows sensor configuration and lsc from pre-loaded (real) meg data

% simulated source coords should be in mm in MEG space

global nuts nutstemp

batch = false
randorient = false
simnoise = true, 
simfilter = false

if(~isfield(nuts,'meg') || ~isfield(nuts.meg,'sensorCoord'))
    errordlg('Load some real MEG data so we can hijack the sensor configuration.');
    return;
end

srate = 1000;   % in Hz
x=-0.750:1/srate:1;  % -500ms:500ms @ 1000 Hz sampling rate

prompt={'source 1:'
        'source 2:'
        'source 3:'
        'SNR (set to ''Inf'' for noiseless signal):'
        'number of trials:'
        'voxel size:' };
name='Simulate Sources';
numlines=[2;2;2;1;1;1];
defaultanswer={strvcat('sin(2*pi*23*(x))','[0 30 40]')
               strvcat('sin(2*pi*23*(x))','[0 -30 40]')
               'none'
               '2'
               '1'
               '5'                                   };

% gammawin = ([zeros(1,550) hanning(100)' zeros(1,400) hanning(120)' zeros(1,331)] + .5)/4;
% gammawin = ([zeros(1,550) hanning(100)' zeros(1,400) hanning(120)' zeros(1,331)]);

taperon = hanning(100)';
taperoff = taperon(51:end);
taperon = taperon(1:50);

% betawin = [ones(1,500) taperoff zeros(1,700) taperon ones(1,201)]+.4;
% alphawin = ([ones(1,600) taperoff zeros(1,650) taperon ones(1,151)]+.6)/2;
% betawin = [ones(1,600) taperoff zeros(1,200) taperon ones(1,601)];
% alphawin = [ones(1,650) taperoff zeros(1,350) taperon ones(1,401)];
           
% defaultanswer={strvcat('sin(2*pi*77*(x+.1)).*gammawin','[0 30 20]')
%                strvcat('sin(2*pi*19*(x-.08)).*betawin','[0 -30 20]')
%                strvcat('sin(2*pi*45*(x-.27)).*alphawin','[-10 0 20]')

simset = 4;
switch(simset)
    case 1 % used for tfsim
            defaultanswer={strvcat('sin(2*pi*77*(x+.1)).*gammawin','[25 30 70]')
                   strvcat('sin(2*pi*19*(x-.08)).*betawin','[-10 40 60]')
                   strvcat('sin(2*pi*45*(x-.27)).*alphawin','[-55 40 50]')
                   '2'
                   '75'
                   '5'                                   };
    case 2 % used for sefa_tf
        defaultanswer={strvcat('sin(2*pi*19*(x+.1)).*(x<0)','[40 50 70]')
                   strvcat('sin(2*pi*77*(x-.27)).*(x>0)','[-55 40 50]')
                   'none'
                   '1'
                   '10'
                   '5'                                   };
    case 3
        gammawin = [zeros(1,650) taperon ones(1,350) taperoff zeros(1,401)];
        betawin = [zeros(1,550) taperon ones(1,300) taperoff zeros(1,200) taperon ones(1,200) taperoff zeros(1,51)];
        defaultanswer={strvcat('sin(2*pi*(77*x+randn)).*gammawin','[40 50 70]')
                   strvcat('sin(2*pi*(25*x+randn)).*betawin','[-55 40 50]')
                   'none'
                   '1'
                   '50'
                   '5'                                   };
    case 4
%         gamma1win = [zeros(1,550) taperon ones(1,350) taperoff zeros(1,501)];
%         gamma2win = [zeros(1,900) taperon ones(1,350) taperoff zeros(1,151)];
%         gamma1win = [zeros(1,750) taperon ones(1,350) taperoff zeros(1,51)];
        gamma1win = [zeros(1,800) taperon ones(1,150) taperoff zeros(1,701)];
        gamma2win = [zeros(1,1100) taperon ones(1,100) taperoff zeros(1,451)];
        betawin = [ones(1,750) taperoff zeros(1,550) taperon ones(1,351)];
        defaultanswer={strvcat('sin(2*pi*(77*x+randn)).*gamma1win','[10 50 60]')
                   strvcat('sin(2*pi*(77*x+randn)).*gamma2win','[15 60 75]')
                   strvcat('sin(2*pi*(19*x+randn)).*betawin','[25 30 100]')
                   '.5'
                   '50'
                   '5'                                   };
end    
 
if(~exist('answer','var'))
    answer=inputdlg(prompt,name,numlines,defaultanswer);
end

if(isempty(answer))
    return; % no input, let's blow this taco stand
end

nuts.params.simparams = answer;


voxels = zeros(0,3);
for i=1:3
    if(size(answer{i},1)==2) % a defined source will have two rows
        % first row contains source; evaluate and add to source list
%         source = [source; eval(answer{i}(1,:))];
        % second row contains coords; add them to voxel list
        voxels = [voxels; str2num(answer{i}(2,:))];
    end
end
realvoxels = nuts.voxels;
nuts.voxels = voxels;
nuts.voxelsize = [5 5 5];

SNR = str2num(answer{4});
numtrials = str2num(answer{5});
voxelsize = str2num(answer{6});

% compute forward solution from specified sources

%% jz's stuff in next 8 lines
% nuts.voxels = voxels;
% % nuts.Lp=nut_compute_lead_field(-1);  % "-1" voxelsize is a hack to get compute_lead_field to cooperate!
% nuts.Lp=nut_compute_lead_field(nuts.voxels);  % now, since this size(input,2)==3, then does simulation hack
% meg_source = (squeeze(nuts.Lp(:,1,:)) * source + squeeze(nuts.Lp(:,2,:)) * source)';  % samples x channels
% meg_noise = (norm(meg_source,'fro')./(SNR*sqrt(prod(size(meg_source)))))*randn(size(meg_source));  % add noise
% meg = meg_source + meg_noise;
% nuts.meg.meg_source=meg_source; % for debugging/verifying, etc.
% nuts.meg.meg_noise=meg_noise; % for debugging/verifying, etc.

tmpvoxels = nuts.voxels;
nuts.voxels = voxels;
nutstemp.Lp=nut_compute_lead_field(voxels);
nuts.voxels = tmpvoxels;
% nuts.voxels = voxels;
% nuts.Lp=nut_compute_lead_field(-1);  % "-1" voxelsize is a hack to get compute_lead_field to cooperate!


L = zeros(size(nutstemp.Lp,1),size(nutstemp.Lp,3));
% note rho and eta are not necessarily consistent with any convention that we use elsewhere!
if(randorient)
%     rho = repmat(rand(1,size(source,1)),size(nuts.Lp,1),1)*2*pi;
    for ii=1:size(nutstemp.Lp,3)
        eta = rand(3,1);
        eta = eta/norm(eta)
        L(:,ii) = nutstemp.Lp(:,:,ii) * eta;
    end
else
%     rho = repmat([.6 1],size(nuts.meg.data,2),1);
%     eta = repmat([.5 .75 .433],size(nuts.meg.data,2),1);
        eta = [1 0 0;0 1 0; 0 0 1]';
    for ii=1:size(nutstemp.Lp,3)
%         eta = [0.5; 0.75; 0.433];
        L(:,ii) = nutstemp.Lp(:,:,ii) * eta(:,ii);
    end
end
nuts.params.eta = eta;

% rho = repmat(.5*ones(1,size(source,1)),size(nuts.Lp,1),1)*2*pi;
% rho = zeros(size(rho));
% meg = (squeeze(nuts.Lp(:,1,:))*source + squeeze(nuts.Lp(:,2,:))*source)';

meg = zeros(size(x,2),size(nutstemp.Lp,1),numtrials);
% zeropoint = dsearchn(x',0);
    
for ii=1:numtrials
    source = zeros(0,length(x));
    for i=1:3
        if(size(answer{i},1)==2) % a defined source will have two rows
            % first row contains source; evaluate and add to source list
            source = [source; eval(answer{i}(1,:))];
        end
    end
    meg(:,:,ii) = (L * source)';  % samples x channels
end

if(simnoise)
    megnoise = zeros(size(meg));
    for ii=1:numtrials
        megnoise(:,:,ii) = (norm(meg(:,:,ii),'fro')./(SNR*sqrt(prod(size(meg(:,:,ii))))))*randn(size(meg(:,:,ii)));
    end
else
    load('megnoise.mat');
    megnoise=megnoise(1:size(meg,1),:,1:size(meg,3));
    % scalefactor dependent on noise dataset:
    % ratio of norm(megnoise(:),'fro') to norm(megnoiseHPF(:),'fro'),
    % where megnoiseHPF is highpass filtered at 2Hz
    scalefactor=1e2;  
    
    noisepower = norm(megnoise(:),'fro');
    megpower = norm(meg(:),'fro');
    SNRorig = megpower/noisepower*scalefactor;
    meg = meg*SNR/SNRorig;  % adjust SNR to desired value
    megpower2 = norm(meg(:),'fro');
    SNR = megpower2/noisepower*scalefactor
end
meg = meg + megnoise;

if(simfilter)
    meg=nut_filter2(meg,'butter','bp',4,2,inf,1000,1);
end

if(~batch)
%     nuts=rmfield(nuts,'voxels'); % remove dummy voxel list used to generate forward solution
    [nuts.Lp,nuts.voxels] = nut_compute_lead_field(voxelsize); % recompute lead field over entire VOI with 5mm voxels
end

%%%% fool the toolbox into thinking this is real.
nuts.meg.data = meg;
nuts.meg.datanoise = megnoise;
nuts.meg.latency = 1000*x';  % convert to milliseconds
nuts.meg.srate = srate;
nuts.meg.filename = 'SimBrain';
if ~isempty(gcbf)  %will be empty during batching
    handles = guihandles(gcbf);
    if isfield(handles,'nut_megfile')
        set(handles.nut_megfile,'String',nuts.meg.filename);
    end
end

nut_enabler;
