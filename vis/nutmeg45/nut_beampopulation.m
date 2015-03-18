 function pop = nut_beampopulation(stats,savefilenames)
% Usages:  pop = nut_beampopulation(stats)
%          nut_beampopulation(stats,savefilenames)
%          nut_beampopulation(stats)
%
%  pop              Output structure containing the population data. If no
%                   output variable is specified, the data is saved to
%                   file.
%  stats            Structure array created with nut_timef_stats GUI.
%  savefilenames    Cell with filenames of population data for each
%                   condition.

% Check input
dosave = (nargout==0);
if dosave
    if nargin>1
        if ~iscell(savefilenames), savefilenames={savefilenames}; end
        if (stats.numtests>1 && length(savefilenames)~=stats.numtests )
            error('You have to specify a savefilename for each condition when using more than 1 condition.')
        end
    else
        if stats.numtests>1
            for k=1:stats.numtests
                [sfile,spath] = uiputfile('*.mat',sprintf('Save population %d as...',k));
                if isequal(sfile,0), return, end
                savefilenames{k} = fullfile(spath,sfile);
            end
        else 
            [sfile,spath] = uiputfile('*.mat','Save population as...');
            if isequal(sfile,0), return, end
            savefilenames{1} = fullfile(spath,sfile);
        end
    end
else
    if numtests>1, error('Can only put 1 condition to the output variable.'), end
end
    
% Get voxelsize and prepare initial parameters
%---------------------------------------------
subj=load(stats.files{stats.subjnr(1,1),stats.condnr(1,1)});  % load first subject to get voxelsize
subj=nut_beam_legacy_compatibility(subj);
voxelsize = subj.voxelsize;

if stats.usespatnormed
    xx = unique(subj.voxels(:,1))';     % Make sure the template grid matches the data voxel grid
    yy = unique(subj.voxels(:,2))';
    zz = unique(subj.voxels(:,3))';
    xx = [ fliplr(xx(1)-voxelsize(1):-voxelsize(1): -90) xx xx(end)+voxelsize(1):voxelsize(1): 90 ];
    yy = [ fliplr(yy(1)-voxelsize(2):-voxelsize(2):-130) yy yy(end)+voxelsize(2):voxelsize(2): 90 ];
    zz = [ fliplr(zz(1)-voxelsize(3):-voxelsize(3): -80) zz zz(end)+voxelsize(3):voxelsize(3):110 ];
    %dimxyz = [length(xx) length(yy) length(zz)];
end

if ( ~stats.docomp && stats.usespatnormed )
    voxels = nut_coordgrid(xx,yy,zz);
    MNItfm = [voxelsize(1) 0 0 xx(1); 0 voxelsize(2) 0 yy(1); 0 0 voxelsize(3) zz(1); 0 0 0 1];
    MNIgrid= uint8(nut_coordtfm(voxels,inv(MNItfm)));  % uint8 = round + liberate memory
else
    voxels=subj.voxels;        % NUT_BEAMCOMP already has all voxels in MNI space
                                    % If processed by NUT_TTEST_FREQBAND, the MNIdeadvoxels are already removed!
end

useroi = isfield(subj,'rois');
if useroi
    roideffile = subj.R.roidef; 
    numroi = length(subj.R.roilabel);
end

numvox = size(voxels,1);
numfiles=size(stats.subjnr,2);

if (isempty(stats.frqband) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.frqband=1:size(subj.bands,1);
end
if (isempty(stats.timepts) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.timepts=1:length(subj.timepts);
end

switch stats.fsel
    case 'one', newfreqbins=length(stats.frqband);
    case 'mean', newfreqbins=1;
end
switch stats.tsel
    case 'one', newtimepts=length(stats.timepts);
    otherwise, newtimepts=1;
end

% Loop across comparisons/tests
%---------------------------------------
for nn=1:stats.numtests            
    if nn>1
        subj=load(stats.files{stats.subjnr(nn,1),stats.condnr(nn,1)});    
        subj=nut_beam_legacy_compatibility(subj);
    end

    % Prepare pop structure
    %-----------------------
    if useroi
        pop.rois  = zeros(stats.numsubj,numroi,length(stats.timepts),newfreqbins);
    else
        pop.s  = zeros(stats.numsubj,numvox,length(stats.timepts),newfreqbins);
    end
    pop.sdims = 'subjects * voxels * timewins * freqbands';
    pop.voxels = voxels;  
    pop.voxelsize = voxelsize;
    pop.timepts = subj.timepts(stats.timepts);           
    pop.timewindow = subj.timewindow(stats.timepts,:);
    switch stats.fsel
        case 'one'
            pop.bands = subj.bands(stats.frqband,:);
        case 'mean'
            pop.bands = [subj.bands(stats.frqband(1),1) subj.bands(stats.frqband(end),2)];
    end
    pop.srate = subj.srate;
    pop.subjnr = stats.subjnr(nn,:);
    pop.condnr = stats.condnr(nn,:);

    %------------------------------
    % Voxels common to all subjects
    %------------------------------
    if useroi, subjgoodroi = zeros(numroi,1); end    
    MNIdeadvoxels = [];
    for ii=1:numfiles     % loop over datasets
        fprintf('Processing subject %d...\n',stats.subjnr(nn,ii))

        %if ii>1     % already loaded otherwise
        subj=load(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)});
        subj=nut_beam_legacy_compatibility(subj);
        if useroi
            subj.s=subj.rois; % for now
            subjgoodroi(subj.R.goodroi)=subjgoodroi(subj.R.goodroi)+1;            
        end 
        if ~isequal(subj.voxelsize,voxelsize), error('NUT_BEAMPOPULATION: All subjects and conditions must have the same voxel size!'), end
        %if useroi && ~strcmp(roideffile,subj.R.roidef), error('NUT_BEAMPOPULATION: All subjects and conditions must have the same ROI definitions.'), end
        if length(subj.s)>1, subj=nut_calc_beamratio(subj); end
        %end

        subj.s{1}=subj.s{1}(:,stats.timepts,stats.frqband);
        if (~stats.frqstats && strcmp(stats.fsel,'mean'))
            subj.s{1}=mean(subj.s{1},3);
        end               
        
        if stats.usespatnormed
            if ~stats.docomp
                voxelgrid = uint8(nut_coordtfm(subj.voxels,inv(MNItfm)));
                [tmp,MNIsubjmap]=ismember(MNIgrid,voxelgrid,'rows');
                clear tmp voxelgrid
                deadvoxels = find(MNIsubjmap==0);
                MNIdeadvoxels = union(MNIdeadvoxels,deadvoxels);
            elseif ~useroi
                if (~stats.frqstats && ~stats.tsstats )
                    subj.s{1}(stats.MNIdeadvoxels{nn},:,:)=[];
                end
                MNIsubjmap=1:size(subj.s{1},1);    % NUT_BEAMCOMP has already arranged subj.s in MNI space!
            end
        else
            MNIsubjmap=1:size(subj.voxels,1);
        end

        %------------------
        % Calculate population
        %------------------     
        if useroi
            pop.rois(ii,subj.R.goodroi,:,:) = subj.s{1};
        else
            %beam.s{1} = zeros(numvox,length(stats.timepts),newfreqbins);
            %for vv=find(MNIsubjmap)'     %1:numvox  % vv num of voxels
            vv=find(MNIsubjmap);
            pop.s(ii,vv,:,:)=subj.s{1}(MNIsubjmap(vv),:,:);
            %beam.s{1}(vv,:,:)=beam.s{1}(vv,:,:) + subj.s{1}(MNIsubjmap(vv),:,:);
            %end
        end
        clear subj MNIsubjmap

    end     % subject loop
    clear voxelgrid deadvoxels

    % remove voxels not present in all subjects or 0 at all time points 
    %------------------------------------------------------------------
    %Q = shiftdim(any(pop.s,1),1);
    %Q = find(any(Q(:,:),2));        
    Q = 1:numvox;                       % Q : voxels non 0 in MNI space
    if stats.usespatnormed
        if stats.docomp 
            Q(stats.MNIdeadvoxels{nn})=[];     
        else    
            Q = setdiff(Q,MNIdeadvoxels)';            
        end
        clear MNIdeadvoxels
    end
    
    pop.voxels=pop.voxels(Q,:); 
    if useroi
        fprintf('Looking for ROI coordinates in common MNI space. This may take a while...\n')
        coreg = struct('mripath','w','norm_mripath','w','meg2mri_tfm',eye(4)); 
        pop.R = fcm_voxel2roi(pop.voxels,coreg,'one',roideffile,find(subjgoodroi<1));
        
        pop.voxels = pop.voxels(pop.R.goodvoxels,:);
        pop.rois = pop.rois(:,pop.R.goodroi,:,:);
        pop.s  = zeros(stats.numsubj,length(pop.R.goodvoxels),length(stats.timepts),newfreqbins);
        for k1=1:stats.numsubj
            for k4=1:newfreqbins
                pop.s(k1,:,:,k4) = pop.R.roi2voxel_tfm * permute(pop.rois(k1,:,:,k4),[2 3 4 1]);
            end
        end
        
    else
        pop.s = pop.s(:,Q,:,:); 
    end
    if dosave
        save(savefilenames{nn},'pop');
        if stats.numtests>1, pop=[]; end
    end
    
    if ~isempty(which('fcm_gui.m'));
        answer=questdlg('Would you like to add this population to a population menu?','NUTMEG question','Patients/Actives','Controls','No thanks','Controls');
        pause(.05) % avoid MATLAB bug which makes the dialog window stay on the screen.
        switch answer
            case 'Patients/Actives'
                answer2=inputdlg('Description of active population:','',1);
                if ~isempty(answer2)
                    load([fileparts(which('fcm.m')) filesep 'params' filesep 'activepopulations.mat']);                  % Add your own control population specs to this file!
                    numexpop = length(actives);
                    actives(numexpop+1).path=savefilenames{nn};
                    actives(numexpop+1).description = answer2{1};
                    save([fileparts(which('fcm.m')) filesep 'params' filesep 'activepopulations.mat'],'actives');
                end
                
            case 'Controls'
                answer2=inputdlg('Description of control population (type, age range, genders):','',1);
                if ~isempty(answer2)
                    load([fileparts(which('fcm.m')) filesep 'params' filesep 'controlpopulations.mat']);                  % Add your own control population specs to this file!
                    numexpop = length(controls);
                    controls(numexpop+1).path=savefilenames{nn};
                    controls(numexpop+1).description = answer2{1};
                    save([fileparts(which('fcm.m')) filesep 'params' filesep 'controlpopulations.mat'],'controls');
                end
        end
    end
               
end     % for nn

fprintf('Done.\n')

