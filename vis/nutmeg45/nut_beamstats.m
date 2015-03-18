 function nut_beamstats(stats)
% NUT_BEAMSTATS: calculates averages and snpm statistics across subjects. 
%
%    Usage:  nut_beamstats(stats)
% stats   Structure array created with nut_timef_stats GUI.
%
% Parts of the code are adapted from the SnPM3 toolbox:  http://www.fil.ion.ucl.ac.uk/spm/snpm
%

% structures used:
% ---------------
% stats:        settings made in nut_timef_stats graphical user interface (GUI).
% SP:           snpm parameters passed to nut_snpm
% beam:         average beamformer data and stat results for display in nut_results_viewer
% subj:         data of single subject.

fprintf('\n\nNUT_BEAMSTATS: Averaging and statistics tool.\n')

% Prepare initial parameters
%---------------------------
if stats.markmode==2
    subj = load(stats.files{1});          % This is pop structure from population file
    if isfield(subj,'pop'), subj = subj.pop; end
else
    subj=loadbeam(stats.files{stats.subjnr(1,1),stats.condnr(1,1)});  % load first subject to get voxelsize
    numfiles=size(stats.subjnr,2);
end

voxelsize = subj.voxelsize;

if stats.usespatnormed
    xx = unique(subj.voxels(:,1))';     % Make sure the template grid matches the data voxel grid
    yy = unique(subj.voxels(:,2))';
    zz = unique(subj.voxels(:,3))';
    xx = [ fliplr(xx(1)-voxelsize(1):-voxelsize(1): -90) xx xx(end)+voxelsize(1):voxelsize(1): 90 ];
    yy = [ fliplr(yy(1)-voxelsize(2):-voxelsize(2):-130) yy yy(end)+voxelsize(2):voxelsize(2): 90 ];
    zz = [ fliplr(zz(1)-voxelsize(3):-voxelsize(3): -80) zz zz(end)+voxelsize(3):voxelsize(3):110 ];
    voxels3d = nut_coordgrid(xx,yy,zz);
else
    voxels3d = subj.voxels;
end
if ( ~stats.docomp && stats.usespatnormed && (stats.markmode<2) )
    MNItfm = [voxelsize(1) 0 0 xx(1); 0 voxelsize(2) 0 yy(1); 0 0 voxelsize(3) zz(1); 0 0 0 1];
    MNIgrid = uint8(nut_coordtfm(voxels3d,inv(MNItfm)));   % uint8 = round + liberate memory
end

if (isempty(stats.frqband) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.frqband=1:size(subj.bands,1);
end
if (isempty(stats.timepts) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.timepts=1:length(subj.timepts);
end

useroi = isfield(subj,'rois');
if useroi
    roideffile = subj.R.roidef; 
    numroi = length(subj.R.roilabel);
end

switch stats.fsel
    case 'one', newfreqbins=length(stats.frqband);
    case 'mean', newfreqbins=1;
end
switch stats.tsel
    case 'one', newtimepts=length(stats.timepts);
    otherwise, newtimepts=1;
end

ftxt=sprintf('_f%s',stats.frqtxt);      % freq info in new filename
if stats.frqstats, ftxt=[ftxt 'stat' stats.whichtail]; end
if strcmp(stats.fsel,'mean'), ftxt=[ftxt 'm']; end
if strcmp(ftxt,'_f'), ftxt=''; end

ttxt=sprintf('_t%s',stats.timetxt);      % time info in new filename
if stats.tsstats, ttxt=[ttxt 'statm']; end
if strcmp(ttxt,'_t'), ttxt=''; end


% Loop across comparisons/tests
%---------------------------------------
for nn=1:stats.numtests          
    
    if length(stats.comps)>1 
        fprintf(['\n\nAnalyzing condition(s) ' num2str(stats.comps(nn,:),' %d') '...\n'])
    end

    if stats.markmode, tmpnr=num2str(stats.comps(nn,:),'%d');
    else, tmpnr=[];
    end
    if stats.doavg, avgfile=['s_beamtf' tmpnr ftxt ttxt '_avg.mat'];
    else, avgfile='';
    end
    if stats.dostat, statfile=['s_beamtf' tmpnr ftxt ttxt '_' strtok(stats.type,'_') '.mat'];
    else, statfile='';
    end
    clear tmpnr
            
    if ( ~exist(avgfile,'file') && ~exist(statfile,'file') )       % skip if file has already been created in earlier session 
        
        % Prepare beam structure
        %-----------------------
        if nn>1
            if stats.markmode<2
                subj=loadbeam(stats.files{stats.subjnr(nn,1),stats.condnr(nn,1)});    
            else
                subj = load(stats.files{nn});
                if isfield(subj,'pop'), subj = subj.pop; end
            end
        end
        beam.timepts = subj.timepts(stats.timepts);           
        beam.timewindow = subj.timewindow(stats.timepts,:);
        switch stats.fsel
        case 'one'
            beam.bands = subj.bands(stats.frqband,:);
        case 'mean'
            beam.bands = [subj.bands(stats.frqband(1),1) subj.bands(stats.frqband(end),2)];
        end
        beam.srate = subj.srate;
        beam.voxelsize=voxelsize;
        if stats.docomp
            beam.voxels = subj.voxels;                     % already have voxels in common space
            if ~(stats.frqstats || stats.tsstats)           % If data was processed by NUT_TTEST_FREQBAND, MNIdeadvoxels are removed already
                beam.voxels(stats.MNIdeadvoxels{nn},:)=[];  % Otherwise, we remove them here
            end
        elseif ( (stats.markmode==2) || ~stats.usespatnormed )
            beam.voxels = subj.voxels;
        else
            beam.voxels = voxels3d; 
        end
        if stats.doavg
            if useroi
                beam.rois{1} = zeros(numroi,length(stats.timepts),newfreqbins);
            else
                beam.s{1} = zeros(size(beam.voxels,1),length(stats.timepts),newfreqbins);
            end
        end
        if stats.usespatnormed || stats.markmode==2
            if strncmp( spm('ver'),'SPM8',4 ) %SPM8 or SPM8b
                spmmri = [fileparts(which('spm')) filesep 'canonical' filesep 'avg152T1.nii'];
            elseif strcmp( spm('ver'),'SPM2' )
                spmmri = [fileparts(which('nutmeg')) filesep 'templates' filesep 'wNormT1neuroax.img'];
            else
                error('This version of SPM is currently not supported.')
            end
            beam.coreg.mripath=spmmri;
            beam.coreg.meg2mri_tfm=eye(4);
            beam.coreg.orientation=1;
        else
            beam.coreg=subj.coreg;
        end
        
        %------------------------------
        % Voxels common to all subjects
        %------------------------------
        if stats.markmode<2
            printflags = [stats.usespatnormed stats.doavg];
            if printflags(1),   fprintf('\nFinding voxels common to all subjects'), end
            if all(printflags), fprintf(' and '), end
            if printflags(2),   fprintf('Calculating average'), end
            clear printflags,   fprintf('...\n')
            
            if useroi, subjgoodroi = zeros(numroi,1); end
            MNIdeadvoxels = [];
            for ii=1:numfiles     % loop over datasets
                fprintf('\tProcessing subject %d',stats.subjnr(nn,ii))
                if ~isempty(strmatch(stats.type,{'anova_snpm','ttest2_snpm'})), fprintf(', condition %d...\n',stats.condnr(ii))
                else fprintf('...\n')
                end

                if ii>1     % already loaded otherwise
                    subj=loadbeam(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)});
                    if ~isequal(subj.voxelsize,voxelsize), error('NUT_BEAMSTATS: All subjects and conditions must have the same voxel size!'), end
                    %if useroi && ~strcmp(roideffile,subj.R.roidef), error('NUT_BEAMSTATS: All subjects and conditions must have the same ROI definitions.'), end                    
                end
                if useroi
                    subj.s=subj.rois; % for now
                    subjgoodroi(subj.R.goodroi)=subjgoodroi(subj.R.goodroi)+1;
                end 
                if length(subj.s)>1, subj=nut_calc_beamratio(subj); end

                subj.s{1}=subj.s{1}(:,stats.timepts,stats.frqband);
                if (~stats.frqstats && strcmp(stats.fsel,'mean'))
                    subj.s{1}=mean(subj.s{1},3);
                end               
                
                if stats.usespatnormed
                    if ~stats.docomp
                        voxelgrid = uint8(nut_coordtfm(subj.voxels,inv(MNItfm)));
                        [tmp,MNIsubjmap]=ismember(MNIgrid,voxelgrid,'rows');
                        clear tmp voxelgrid
                        save(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)},'-append','MNIsubjmap');
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
                % Calculate average
                %------------------     
                if stats.doavg
                    if useroi
                        beam.rois{1}(subj.R.goodroi,:,:) = beam.rois{1}(subj.R.goodroi,:,:) + subj.s{1};
                    else
                        %for jj=1:length(subj.s)  %jj num of orientations
                            %%%%%%%%% needs lots of RAM for this method
                        good = (MNIsubjmap>0);
                        beam.s{1}(good,:,:) = beam.s{1}(good,:,:) + subj.s{1}(MNIsubjmap(good),:,:);
                            %%%%%%%%% this is slower but needs far less RAM
                            % for vv=find(MNIsubjmap)     %1:numvox  % vv num of voxels
                            %    beam.s{jj}(vv,:,:)=beam.s{jj}(vv,:,:) + subj.s{jj}(MNIsubjmap(vv),:,:);
                            % end
                        %end
                    end
                end
                %         pop.s  = zeros(stats.numsubj,length(Q),length(stats.timepts),newfreqbins);
                %         for k1=1:stats.numsubj
                %             for k4=1:newfreqbins
                %                 pop.s(k1,:,:,k4) = R.roi2voxel_tfm * permute(pop.rois(k1,:,:,k4),[2 3 4 1]);
                %             end
                %         end
                clear subj MNIsubjmap good

            end     % subject loop
            clear voxelgrid deadvoxels
            if stats.doavg
                if useroi
                    beam.rois{1} = beam.rois{1}./stats.numsubj;  
                else
                    beam.s{1} = beam.s{1}./stats.numsubj;  
                end
            end
                        
        else        % stats.markmode==2
            if useroi
                beam.rois{1} = shiftdim(mean(subj.rois(:,:,stats.timepts,stats.frqband),1),1);
                beam.R = subj.R;
            else
                beam.s{1} = shiftdim(mean(subj.s(:,:,stats.timepts,stats.frqband),1),1); 
            end
            clear subj
        end 
        
        % remove voxels not present in all subjects or 0 at all time points 
        %------------------------------------------------------------------
        Q = 1:size(voxels3d,1);
        if stats.docomp 
            Q(stats.MNIdeadvoxels{nn})=[];          % Q: voxels defined in all subjects
        else
            if stats.markmode==2
                MNIdeadvoxels = find(~ismember(voxels3d,beam.voxels,'rows'));
            end
            if stats.usespatnormed
                Q(MNIdeadvoxels)=[];
            end
            if ( (stats.markmode<2) && stats.usespatnormed )
                beam.voxels = voxels3d(Q,:);
                if stats.doavg && ~useroi
                    beam.s{1}=beam.s{1}(Q,:,:);
                end
            end
        end
        beam.MNIgoodvoxels=Q;
        
        if useroi
            if ~isfield(beam,'R')
                fprintf('\tLooking for ROI coordinates in common MNI space. This may take a while...\n')
                coreg = struct('mripath','w','norm_mripath','w','meg2mri_tfm',eye(4)); 
                beam.R = fcm_voxel2roi(beam.voxels,coreg,'one',roideffile,find(subjgoodroi<1));

                beam.voxels = beam.voxels(beam.R.goodvoxels,:);
                if stats.doavg, beam.rois{1} = beam.rois{1}(beam.R.goodroi,:,:); end
            end
            if stats.doavg
                beam.rois{1}(isnan(beam.rois{1}))=0;
                beam.s{1}  = zeros(length(beam.R.goodvoxels),length(stats.timepts),newfreqbins);
                for k3=1:newfreqbins
                    beam.s{1}(:,:,k3) = beam.R.roi2voxel_tfm * beam.rois{1}(:,:,k3);
                end
            end
        end

        % Save data calculated so far
        %----------------------------
        if stats.doavg
            fprintf('\tSaving average as %s...\n',avgfile)
            save(avgfile,'-struct','beam');
        elseif stats.dostat
            fprintf('\tSaving voxel structure as %s...\n',statfile)
            save(statfile,'-struct','beam');
        end

    elseif stats.dostat   % if a file already exists and stats requested, load average/voxel structure calculated in former session
        if stats.doavg
            fprintf('\nUsing average file %s generated in earlier session.\n',avgfile)
            beam=loadbeam(avgfile);
        else
            fprintf('\nUsing file %s generated in earlier session.\n',statfile)
            beam=loadbeam(statfile);            
        end
        Q=beam.MNIgoodvoxels;
    end
    
    
    %---------------------------------------
    % Stats part: SnPM 
    %---------------------------------------
    if stats.dostat  % if stats requested
        
        if isfield(beam,'snpm')
            if ( (isfield(beam.snpm,'p_corr_pos') && ~isempty(strmatch(stats.whichtail,{'pos','both'}))) || (isfield(beam.snpm,'p_corr_neg') && ~isempty(strmatch(stats.whichtail,{'neg','both'}))) )
                answer=questdlg('WARNING: Results of SnPM tests with the current hypothesis are already present in your current s_beam*.mat file. If you continue, these results will be overwritten and lost. Do you wish to continue?', ...
                    'NUTMEG Question','Yes, go ahead','Oh no, stop please','Yes, go ahead');
                pause(.1)   % This avoids stupid Matlab bug that freezes all display till calcs are done.
                if strcmp(answer,'Oh no, stop please'), error('NUT_BEAMSTATS: User abort.'), end
            end
        end
        %if stats.docomp, MNI0voxels=beam.MNI0voxels; end
        %beam = struct('voxels',beam.voxels,'voxelsize',beam.voxelsize,'coreg',beam.coreg);      % save memory
        
        SP = make_permmatrix(stats);
        
        if useroi
            SP.dosmoothvar=false;
            numgvox=numroi; 
        else
            numgvox=length(Q); 
        end
        
        nP   = ones(numgvox,newtimepts,newfreqbins);        % initialize
        T    = zeros(numgvox,newtimepts,newfreqbins);
        MaxT = zeros(SP.numperm,(SP.STAT~='F')+1,newfreqbins,newtimepts);       % permutations*[max min]*frqbands*timepts
        MaxCluster = zeros(SP.numperm,(SP.STAT~='F')+1,newfreqbins,newtimepts);
        
        % Loop through time windows and frequency bands
        %----------------------------------------------  
        for tt=1:newtimepts            % loop through time windows if requested.
            if newtimepts>1, fprintf('\n\tPROCESSING TIME WINDOW %d\n',stats.timepts(tt)), end
            for ff=1:newfreqbins           % loop through frequency bands if requested.
                if newfreqbins>1, fprintf('\n\tPROCESSING FREQENCY BAND %d\n',stats.frqband(ff)), end

                fprintf('\tLoading data...\n')
                switch stats.fsel
                    case 'one'
                        fprintf('\t\tFrequency band setting: current.\n')
                    case 'mean'
                        fprintf('\t\tFrequency band setting: average ')
                        if stats.frqstats, fprintf('of statistically significant.\n')
                        else, fprintf('of all selected.\n')
                        end
                    otherwise
                        error('NUT_BEAMSTATS: Your stats settings were messed up.')
                end
                switch stats.tsel
                    case 'meanstat'
                        fprintf('\t\tTime window setting: average of statistically significant.\n')
                    case 'one'
                        fprintf('\t\tTime window setting: current.\n')
                    case 'sum'
                        fprintf('\t\tTime window setting: area under curve of all selected.\n')
                    otherwise
                        error('NUT_BEAMSTATS: Your stats settings were messed up.')
                end
                
                if stats.markmode<2
                    X = zeros(numfiles,numgvox);       % numfiles = num_subjects*num_conditions in case of ANOVA
                    for ii=1:numfiles    % loop over datasets
                        subj = loadbeam(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)});
                        if useroi, subj.s=subj.rois; end % for now
                        if length(subj.s)>2, subj=nut_calc_beamratio(subj); end
                        if stats.docomp && ~useroi
                            if ( ~stats.frqstats && ~stats.tsstats )
                                subj.s{1}(stats.MNIdeadvoxels{nn},:,:)=[];
                            end
                            %subj.s{1}(MNI0voxels,:,:)=[];
                        end
                        switch stats.fsel
                            case 'one'
                                subj.s{1}=subj.s{1}(:,stats.timepts,stats.frqband(ff));
                            case 'mean'
                                subj.s{1}=mean(subj.s{1}(:,stats.timepts,stats.frqband),3);
                        end
                        switch stats.tsel
                            case 'one'
                                ts=subj.s{1}(:,tt);
                            case 'meanstat'
                                ts = zeros(size(subj.s{1},1),1);
                                for vv=1:size(subj.s{1},1)
                                    if ~isempty(subj.goodtimepts{vv})
                                        ts(vv) = mean( subj.s{1}(vv,subj.goodtimepts{vv}) ,2 );
                                    end
                                end
                            case 'sum'
                                %ts=sum(subj.s{1},2);
                                ts=trapz(subj.timepts./1000,subj.s{1},2);   % Calculate area under curve (AUC)
                        end
                        
                        if useroi
                            X(ii,subj.R.goodroi)=ts;
                        elseif ( stats.docomp || ~stats.usespatnormed )
                            X(ii,:)=ts;
                        else
                            for vv=1:numgvox  
                                X(ii,vv) = ts(subj.MNIsubjmap(Q(vv)));
                            end
                        end
                        clear ts subj
                    end  % subject loop
                    
                else     %if stats.markmode==2
                    pop = load(stats.files{nn});
                    if isfield(pop,'pop'), pop=pop.pop; end
                    if useroi, pop.s= pop.rois; end % for now
                    if isfield(stats,'subj2remove')
                        pop.s(stats.subj2remove,:,:,:)=[]; 
                    else
                        pop.s=pop.s(ismember(pop.subjnr,stats.subjnr(nn,:)),:,:,:);
                    end
                    switch stats.fsel
                        case 'one'
                            pop.s=pop.s(:,:,stats.timepts,stats.frqband(ff));
                        case 'mean'
                            pop.s=mean(pop.s(:,:,stats.timepts,stats.frqband),4);
                    end
                    switch stats.tsel
                        case 'one'
                            X=pop.s(:,:,tt,:);
                        case 'meanstat'
                            error('The option "Use mean of time points which are significant across frequency bands" is not implemented for population files. Use a pointer file instead.')
                        case 'sum'
                            X=trapz(pop.timepts./1000,pop.s,3);   % Calculate area under curve (AUC)
                    end
                    clear pop
                end
                
                % SnPM main part (permutations)
                %------------------------------
                [T(:,tt,ff),nP(:,tt,ff),MaxT(:,:,ff,tt),MaxCluster(:,:,ff,tt),cth1] = nut_snpm(X,Q,voxels3d,voxelsize,SP);

            end     % freqband loop
        end     % timepts loop
        if all(isnan(MaxCluster(:))), MaxCluster=[]; end
        if newfreqbins>1 || newtimepts>1, fprintf('\n'), end
       
        % p values
        %------------
        fprintf('\tCalculating p-values...\n')
        nP = nP/stats.numperm; 

        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            nP_pos=nP;
            nP_pos(isnan(T))=1;
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT~='F'
            nP_neg = 1+1/stats.numperm-nP;
            nP_neg((T==0) | isnan(T)) = 1;
        end
        clear nP

        if newfreqbins>1 
            if stats.corr4multfreq
                fprintf('\tCorrecting for testing multiple frequency bands.\n')
                MaxT = repmat(max(MaxT,[],3),[1 1 newfreqbins 1]);          % correct for multiple frequency bands
                MaxCluster = repmat(max(MaxCluster,[],3),[1 1 newfreqbins 1]);
            else
                fprintf('\tDoes not correct for testing multiple frequency bands.\n')
            end
        end
        if newtimepts>1 
            if stats.corr4multtime
                fprintf('\tCorrecting for testing multiple time windows.\n')
                MaxT = repmat(max(MaxT,[],4),[1 1 1 newtimepts]);          % correct for multiple time windows
                MaxCluster = repmat(max(MaxCluster,[],4),[1 1 1 newtimepts]);
            else
                fprintf('\tDoes not correct for testing multiple time windows.\n')
            end
        end

        if SP.hPerms
            MaxT = cat( 1, MaxT, flipdim(flipdim(MaxT,2),1) );
            MaxCluster = cat(1, MaxCluster, flipdim(flipdim(MaxCluster,2),1) );
        end
        
        tol = 1e-4;	% Tolerance for comparing real numbers
        iFWE=ceil(.95*stats.numperm);
        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            % Voxel correction
            cP_pos = zeros(numgvox,newtimepts,newfreqbins);    
            for pp=1:stats.numperm    %MaxT_pos'
                %-FEW-corrected p is proportion of randomisation greater or
                % equal to statistic.
                cP_pos = cP_pos + ( repmat(permute(MaxT(pp,1,:,:),[1 4 3 2]),[numgvox 1 1]) > T-tol );      % correct for multiple voxels
            end
            cP_pos = cP_pos / stats.numperm;  
            cP_pos(isnan(T)) = 1;

            % Cluster correction
            if ~isempty(MaxCluster)
                C_pos = zeros( newtimepts - ( isfield(stats,'corr4multtime') && stats.corr4multtime )*(newtimepts-1) , newfreqbins - ( isfield(stats,'corr4multfreq') && stats.corr4multfreq )*(newfreqbins-1));
                for ff=[1 : newfreqbins-( isfield(stats,'corr4multfreq') && stats.corr4multfreq )*(newfreqbins-1)]
                    for tt=[1 : newtimepts-( isfield(stats,'corr4multtime') && stats.corr4multtime )*(newtimepts-1)]
                        MCS = sort(MaxCluster(:,1,ff,tt));
                        C_pos(tt,ff) = MCS(iFWE);
                    end
                end
                sC = size(C_pos);
                if ~isequal(sC,[newtimepts newfreqbins]) && ~isscalar(C_pos)
                    C_pos = repmat(C_pos,[newtimepts newfreqbins]./sC);
                end
                cPclust_pos = nut_clusterstats(T,nP_pos,beam,C_pos,cth1);
            end
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT~='F'
            cP_neg = zeros(numgvox,newtimepts,newfreqbins);
            for pp = 1:stats.numperm    %MaxT_neg'
                cP_neg = cP_neg + ( repmat(permute(MaxT(pp,2,:,:),[1 4 3 2]),[numgvox 1 1]) > -T-tol );
            end
            cP_neg = cP_neg / stats.numperm;  
            cP_neg(isnan(T)) = 1;
            
            % Cluster correction
            if ~isempty(MaxCluster)
                C_neg = zeros( newtimepts - ( isfield(stats,'corr4multtime') && stats.corr4multtime )*(newtimepts-1) , newfreqbins - ( isfield(stats,'corr4multfreq') && stats.corr4multfreq )*(newfreqbins-1));
                for ff=[1 : newfreqbins-( isfield(stats,'corr4multfreq') && stats.corr4multfreq )*(newfreqbins-1)]
                    for tt=[1 : newtimepts-( isfield(stats,'corr4multtime') && stats.corr4multtime )*(newtimepts-1)]
                        MCS = sort(MaxCluster(:,2,ff,tt));
                        C_neg(tt,ff) = MCS(iFWE);
                    end
                end
                sC = size(C_neg);
                if ~isequal(sC,[newtimepts newfreqbins]) && ~isscalar(C_neg)
                    C_neg = repmat(C_neg,[newtimepts newfreqbins]./sC);
                end
                cPclust_neg = nut_clusterstats(-T,nP_neg,beam,C_neg,cth1);
            end
        end
        clear MaxT MaxCluster

        if useroi
            fprintf('\tProjecting ROI data to voxels...\n')
            RT = T(beam.R.goodroi,:,:); 
            RT(isnan(RT))=0;
            T  = zeros(length(beam.R.goodvoxels),length(stats.timepts),newfreqbins); 
            for k3=1:newfreqbins
                T(:,:,k3) = beam.R.roi2voxel_tfm * RT(:,:,k3);
            end
            todolist = {'nP_pos' 'nP_neg' 'cP_pos' 'cP_neg'};
            todolist = todolist([exist('nP_pos','var') exist('nP_neg','var') exist('cP_pos','var') exist('cP_neg','var')]==1);
            for ll=1:length(todolist)
                eval(['R' todolist{ll} '=' todolist{ll} '(beam.R.goodroi,:,:);']);
                eval([todolist{ll} '= zeros(length(beam.R.goodvoxels),length(stats.timepts),newfreqbins);']);
                for k3=1:newfreqbins
                    eval([todolist{ll} '(:,:,k3)= beam.R.roi2voxel_tfm * R' todolist{ll} '(:,:,k3);']);
                end
            end
            clear todolist
        end
        
        % Save stats in beam structure
        %------------------------------------
%         if stats.doavg
%             beam=loadbeam(avgfile); 
%         else
%             beam=loadbeam(statfile);
%         end
        
        beam.snpm.T=T;
        if useroi, beam.snpm.roi_T=RT; end
        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            beam.snpm.p_uncorr_pos=nP_pos;
            beam.snpm.p_corr_pos=cP_pos;
            if useroi
                beam.snpm.roi_p_uncorr_pos=RnP_pos;
                beam.snpm.roi_p_corr_pos=RcP_pos;
            end
            if exist('cPclust_pos','var'), beam.snpm.p_cluster_corr_pos=cPclust_pos; end
            if exist('C_pos','var'), beam.snpm.clusterthres_pos = C_pos; end
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT~='F'
            beam.snpm.p_uncorr_neg=nP_neg;
            beam.snpm.p_corr_neg=cP_neg;
            if useroi
                beam.snpm.roi_p_uncorr_neg=RnP_neg;
                beam.snpm.roi_p_corr_neg=RcP_neg;
            end            
            if exist('cPclust_neg','var'), beam.snpm.p_cluster_corr_neg=cPclust_neg; end
            if exist('C_neg','var'), beam.snpm.clusterthres_neg = C_neg; end
        end   
        % add settings structure
        beam.snpm.settings.model=strtok(stats.type,'_');     
        beam.snpm.settings.subj.number=stats.numsubj;
        beam.snpm.settings.subj.selected=unique(stats.subjsel);
        beam.snpm.settings.freq.analysis=stats.fsel;
        beam.snpm.settings.freq.stats=stats.frqstats;
        if isfield(stats,'corr4multfreq'), beam.snpm.settings.freq.correct=stats.corr4multfreq; end
        beam.snpm.settings.time.analysis=stats.tsel;
        if isfield(stats,'corr4multime'), beam.snpm.settings.time.correct=stats.corr4multtime; end
        
        if stats.doavg
            fprintf('\tAdding statistics to average file %s...\n',avgfile)       
            save(avgfile,'-struct','beam');           % resave average file with stats data, no masking of nonsignificant voxels
        end

        fprintf('\tSaving statistics as %s...\n',statfile)
        beam.s{1}=beam.snpm.T;
        if useroi, beam.rois{1}=beam.snpm.roi_T; end
        beam.sinfo={SP.STAT};
        if ~strcmp(stats.tsel,'one')
            beam.timepts = [beam.timepts(1);beam.timepts(end)];                 % Adjust time info if not stats at each time window calculated
            beam.timewindow = [beam.timewindow(1,1) beam.timewindow(end,2)];
        end
        save(statfile,'-struct','beam');

        clear beam Q T nP* cP* avgfile statfile
    
        
    else
        clear beam Q
        
    end     % if stats.dostat
end     % numtests loop

fprintf('Done.\n')


%====================================================================================
function SP = make_permmatrix(stats)

switch stats.type 

% t-test for 1 sample
%---------------------
case {'snpm','ttest1_snpm'}     

    fprintf('\nPerforming SnPM t-tests for one sample across subjects...\n')      

    SP.sHCform  = 'spm_DesMtx(SP.PiCond(perm,:),''C'',''Mean'')';
    SP.CONT     = 1;
    SP.STAT     = 'T';            
    iCond = ones(1,stats.numsubj);

    % create permutation matrix
    %--------------------------
    fprintf('\tCreating permutation matrix...\n')

    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.hPerms = ~stats.isaproxperm;     %-Only do half the work, if possible
    if stats.numsubj<=12 || ~stats.isaproxperm
        SP.PiCond=[];
        for i=0:stats.numsubj-1
            SP.PiCond=[ones(2^i,1),SP.PiCond;-ones(2^i,1),SP.PiCond];
        end
        if ~stats.isaproxperm    
            SP.PiCond=SP.PiCond(SP.PiCond(:,1)==1,:);
            SP.numperm=stats.numperm/2;            % note: SP.numperm is number of permutations actually calculated, stats.numperm the total number
        else
            SP.numperm=stats.numperm;
            tmp=randperm(size(SP.PiCond,1));
            SP.PiCond=SP.PiCond(tmp(1:SP.numperm),:);
            % Note we may have missed iCond!  We catch this below.	
        end
    else
        SP.numperm=stats.numperm;
        d       = SP.numperm-1;
        tmp     = pow2(0:stats.numsubj-1)*ones(stats.numsubj,1);  % Include correctly labeled iCond

        while (d>0)
          tmp = union(tmp,floor(rand(1,d)*2^stats.numsubj));
          tmp(tmp==2^stats.numsubj) = [];  % This will almost never happen
          d   = SP.numperm-length(tmp);
        end

        % randomize tmp before it is used to get SP.PiCond
        rand_tmp=randperm(length(tmp));
        tmp=tmp(rand_tmp);
        SP.PiCond = 2*rem(floor(tmp(:)*pow2(-(stats.numsubj-1):0)),2)-1;
    end
    clear tmp d rand_tmp i

    %-Find (maybe) iCond in PiCond, move iCond to 1st; negate if neccesary
    %-----------------------------------------------------------------------
    perm = find(all((meshgrid(iCond,1:size(SP.PiCond,1))==SP.PiCond)'));
    if SP.hPerms
        perm=[perm,-find(all((meshgrid(iCond,1:size(SP.PiCond,1))==-SP.PiCond)'))];
    end
    if length(perm)==1
        if (perm<0), SP.PiCond=-SP.PiCond; perm=-perm; end
        %-Actual labelling must be at top of PiCond
        if (perm~=1)
            SP.PiCond(perm,:)=[];
            SP.PiCond=[iCond;SP.PiCond];
        end
        if ~stats.isaproxperm 
            %-Randomise order of PiConds, unless already randomized
            % Allows interim analysis	
            SP.PiCond=[SP.PiCond(1,:);SP.PiCond(randperm(size(SP.PiCond,1)-1)+1,:)];
        end	
    elseif isempty(perm) && (stats.numsubj<=12) && stats.isaproxperm
        % Special case where we missed iCond; order of perms is random 
        % so can we can just replace first perm.
        SP.PiCond(1,:) = iCond;
    else    
        error(['Bad PiCond (' num2str(perm) ')'])
    end    

% unpaired t-test
%-------------------------------------------------
case 'ttest2_snpm'      

    fprintf('\nPerforming SnPM unpaired t-tests across subjects...\n')      
    
    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.sHCform = 'spm_DesMtx(SP.PiCond(perm,:),''-'',''Cond'')';
    SP.CONT = [-1 1];        
    SP.STAT     = 'T';

    % create permutation matrix
    %--------------------------
    fprintf('\tCreating permutation matrix...\n')  

    rand('seed',sum(100*clock));        % Initialise random number generator
    iCond = stats.condnr-min(stats.condnr); iCond = -iCond/max([1,iCond])*2+1;
    nFlip=length(find(iCond==-1));    
    SP.hPerms = (~stats.isaproxperm && (nFlip == stats.numsubj/2));   % balanced group numbers

    if stats.numsubj<=12 || ~stats.isaproxperm        % exact method
        SP.PiCond=[];
        for i=0:stats.numsubj-1     %-Generate all labellings of stats.numsubj scans as +/- 1
            SP.PiCond=[ones(2^i,1),SP.PiCond;-ones(2^i,1),SP.PiCond];
        end    
        SP.PiCond=SP.PiCond(sum(SP.PiCond'==-1)==nFlip,:);  %-Trim to labellings with correct group numbers
        SP.numperm=stats.numperm;
        if SP.hPerms  % do half the work if possible
            SP.PiCond=SP.PiCond(SP.PiCond(:,1)==1,:);
            SP.numperm=stats.numperm/2;
        elseif stats.isaproxperm        % pick random supsample of perms                                   
            tmp=randperm(size(SP.PiCond,1));
            PiCond=PiCond(tmp(1:stats.numperm),:);
            % Note we may have missed iCond!  We catch this below.	
        end

    else                                          % random method
        SP.numperm=stats.numperm;
        SP.PiCond = zeros(stats.numperm,stats.numsubj);    % Allocate final result
        SP.PiCond(1,:) = iCond;                    % Fill first row 
        for i=2:stats.numperm                           % Fill subsequent rows, checking that we're not repeating
            tmp=SP.PiCond(i-1,randperm(stats.numsubj));
            while any(all(SP.PiCond(1:(i-1),:)'==meshgrid(tmp,1:(i-1))'))
                tmp=SP.PiCond(i-1,randperm(stats.numsubj));
            end
            SP.PiCond(i,:)=tmp;
        end      
    end        
    if ~all(all(SP.PiCond*ones(stats.numsubj,1)==stats.numsubj-2*nFlip))    %Check PiConds sum to nGrp1-nGrp2
        error('Invalid PiCond computed!')
    end

    %-Find (maybe) iCond in SP.PiCond, move iCond to 1st; negate if neccesary
    perm = find(all((meshgrid(iCond,1:size(SP.PiCond,1))==SP.PiCond)'));
    if SP.hPerms
        perm=[perm,-find(all((meshgrid(iCond,1:size(SP.PiCond,1))==-SP.PiCond)'))];
    end
    if length(perm)==1
        if (perm<0), SP.PiCond=-SP.PiCond; perm=-perm; end
        %-Actual labelling must be at top of PiCond
        if (perm~=1)
            SP.PiCond(perm,:)=[];
            SP.PiCond=[iCond;SP.PiCond];
        end
        if ~stats.isaproxperm
            %-Randomise order of PiConds, unless already randomized
            % Allows interim analysis	
            SP.PiCond=[SP.PiCond(1,:);SP.PiCond(randperm(size(SP.PiCond,1)-1)+1,:)];
        end
    else    
        error(['Bad PiCond (' num2str(perm) ')'])
    end    
    clear tmp i perm

% ANOVA
%------------------------------
case 'anova_snpm'

    fprintf('\nPerforming SnPM ANOVA tests across subjects...\n')      

    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.iSubj    = stats.subjnr;
    SP.iRepl    = stats.condnr;
    SP.sHCform  = 'diag(SP.PiCond(perm,:)*spm_DesMtx(SP.iSubj)'') * spm_DesMtx(SP.iRepl,''-'',''Scan'')';
    SP.CONT     = eye(length(stats.comps));
    SP.STAT     = 'F';

    % create permutation matrix
    %--------------------------
    fprintf('\tCreating permutation matrix...\n') 

    iCond = ones(1,stats.numsubj);
    % Set SP.hPerms=0. The reason is this:
    % the permutations with all +1's or all -1's will give the same F.
    % So we just want to count half of all possible permutations.
    % Another way to think about it is to always keep first subject as +1.
    SP.hPerms=false;
    SP.numperm=stats.numperm;

    if stats.numsubj<=12 || ~stats.isaproxperm          % exact method

        %-Generate all labellings of stats.numsubj scans as +/- 1
        SP.PiCond=[];
        for i=0:stats.numsubj-2
            SP.PiCond=[ones(2^i,1),SP.PiCond;-ones(2^i,1),SP.PiCond];
        end

        a = ones(size(SP.PiCond,1),1);
        SP.PiCond =[a,SP.PiCond];

        if stats.isaproxperm                 % pick random supsample of perms
            tmp=randperm(size(SP.PiCond,1));
            SP.PiCond=SP.PiCond(tmp(1:nPiCond),:);
                % Note we may have missed iCond!  We catch this below.	
        end	

    else                                          % random method
        d       = stats.numperm-1;
        tmp     = pow2(0:stats.numsubj-2)*iCond(1:(stats.numsubj-1))';  % Include correctly labeled iCond

        while (d>0)
           tmp = union(tmp,floor(rand(1,d)*2^(stats.numsubj-1)));
           tmp(tmp==2^(stats.numsubj-1)) = [];  % This will almost never happen
           d   = stats.numperm-length(tmp);
        end

        % randomize tmp before it is used to get SP.PiCond
        rand_tmp=randperm(length(tmp));
        tmp=tmp(rand_tmp);

        SP.PiCond = 2*rem(floor(tmp(:)*pow2(-(stats.numsubj-1-1):0)),2)-1;

        a = ones(size(SP.PiCond,1),1);
        SP.PiCond =[a,SP.PiCond]; 
    end

    %-Find (maybe) iCond in SP.PiCond, move iCond to 1st; negate if neccesary
    %-----------------------------------------------------------------------
    perm = find(all((meshgrid(iCond,1:size(SP.PiCond,1))==SP.PiCond)'));
    if length(perm)==1
        if (perm<0), SP.PiCond=-SP.PiCond; perm=-perm; end
        %-Actual labelling must be at top of SP.PiCond
        if (perm~=1)
            SP.PiCond(perm,:)=[];
            SP.PiCond=[iCond;SP.PiCond];
        end
        if ~stats.isaproxperm    
            %-Randomise order of PiConds, unless already randomized
            % Allows interim analysis	
            SP.PiCond=[SP.PiCond(1,:);SP.PiCond(randperm(size(SP.PiCond,1)-1)+1,:)];
        end	
    elseif length(perm)==0 && (stats.numsubj<=12) && stats.isaproxperm
        % Special case where we missed iCond; order of perms is random 
        % so can we can just replace first perm.
        SP.PiCond(1,:) = iCond;
        perm = 1;
    else    
        error(['Bad PiCond (' num2str(perm) ')'])
    end                      
    clear tmp a d rand_tmp perm
    
% Across subject correlation
%------------------------------
case 'corr_snpm'
    fprintf('\nPerforming SnPM correlations across subjects...\n')      

    SP.sHCform  = ['spm_DesMtx( SP.PiCond(perm,:) .* SP.CovInt.'',''C'',''CovInt(centred)'')'];
    SP.CONT     = [1];        
    answer = questdlg('Which stats do you want to calculate from your correlations?','NUTMEG question','Variance-smoothed pseudo-T','Non-variance-smoothed Pearson R','Non-variance-smoothed Pearson R');
    pause(.1)   % This avoids stupid Matlab bug that freezes display till calcs are done.
    if strcmp(answer,'Variance-smoothed pseudo-T')
        SP.STAT = 'T';
    else
        SP.STAT = 'R';
    end
    SP.CovInt   = (stats.behavdata-mean(stats.behavdata)) ./ norm(stats.behavdata-mean(stats.behavdata));
    
    iCond = ones(1,stats.numsubj);

    % create permutation matrix
    %--------------------------
    fprintf('\tCreating permutation matrix...\n')

    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.hPerms = ~stats.isaproxperm;     %-Only do half the work, if possible
    if stats.numsubj<=12 || ~stats.isaproxperm
        SP.PiCond=[];
        for i=0:stats.numsubj-1
            SP.PiCond=[ones(2^i,1),SP.PiCond;-ones(2^i,1),SP.PiCond];
        end
        if ~stats.isaproxperm    
            SP.PiCond=SP.PiCond(SP.PiCond(:,1)==1,:);
            SP.numperm=stats.numperm/2;            % note: SP.numperm is number of permutations actually calculated, stats.numperm the total number
        else
            SP.numperm=stats.numperm;
            tmp=randperm(size(SP.PiCond,1));
            SP.PiCond=SP.PiCond(tmp(1:SP.numperm),:);
            % Note we may have missed iCond!  We catch this below.	
        end
    else
        SP.numperm=stats.numperm;
        d       = SP.numperm-1;
        tmp     = pow2(0:stats.numsubj-1)*ones(stats.numsubj,1);  % Include correctly labeled iCond

        while (d>0)
          tmp = union(tmp,floor(rand(1,d)*2^stats.numsubj));
          tmp(tmp==2^stats.numsubj) = [];  % This will almost never happen
          d   = SP.numperm-length(tmp);
        end

        % randomize tmp before it is used to get SP.PiCond
        rand_tmp=randperm(length(tmp));
        tmp=tmp(rand_tmp);
        SP.PiCond = 2*rem(floor(tmp(:)*pow2(-(stats.numsubj-1):0)),2)-1;
    end
    clear tmp d rand_tmp i

    %-Find (maybe) iCond in PiCond, move iCond to 1st; negate if neccesary
    %-----------------------------------------------------------------------
    perm = find(all((meshgrid(iCond,1:size(SP.PiCond,1))==SP.PiCond)'));
    if SP.hPerms
        perm=[perm,-find(all((meshgrid(iCond,1:size(SP.PiCond,1))==-SP.PiCond)'))];
    end
    if length(perm)==1
        if (perm<0), SP.PiCond=-SP.PiCond; perm=-perm; end
        %-Actual labelling must be at top of PiCond
        if (perm~=1)
            SP.PiCond(perm,:)=[];
            SP.PiCond=[iCond;SP.PiCond];
        end
        if ~stats.isaproxperm 
            %-Randomise order of PiConds, unless already randomized
            % Allows interim analysis	
            SP.PiCond=[SP.PiCond(1,:);SP.PiCond(randperm(size(SP.PiCond,1)-1)+1,:)];
        end	
    elseif isempty(perm) && (stats.numsubj<=12) && stats.isaproxperm
        % Special case where we missed iCond; order of perms is random 
        % so can we can just replace first perm.
        SP.PiCond(1,:) = iCond;
    else    
        error(['Bad PiCond (' num2str(perm) ')'])
    end      
    
case 'corrorig_snpm'    
    fprintf('\nPerforming SnPM correlations across subjects...\n')      

    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.sHCform  = ['spm_DesMtx( (SP.CovInt(SP.PiCond(perm,:))-mean(SP.CovInt)) ./ norm(SP.CovInt(SP.PiCond(perm,:))-mean(SP.CovInt)) ,',...
                    '''C'',''CovInt(centred)'')'];
    SP.CONT     = [1];        
    SP.STAT     = 'R';
    SP.CovInt   = stats.behavdata;
        
    % create permutation matrix
    %--------------------------
    fprintf('\tCreating permutation matrix...\n')     

    SP.hPerms=false;
    SP.numperm=stats.numperm;
    if stats.isaproxperm
        %-Approximate test :
        % Build up random subset of all (within nSubj) permutations
        %===============================================================
        rand('seed',sum(100*clock))	%-Initialise random number generator
        SP.PiCond      = zeros(stats.numperm,stats.numsubj);
        SP.PiCond(1,:) = 1+rem([0:stats.numsubj-1],stats.numsubj);
        for i = 2:stats.numperm
            %-Generate a new random permutation - see randperm
            [null,p] = sort(rand(stats.numsubj,1)); p = p(:)';
            %-Check it's not already in SP.PiCond
            while any(all((meshgrid(p,1:i-1)==SP.PiCond(1:i-1,:))'))
                [null,p] = sort(rand(stats.numsubj,1)); p = p(:)';
            end
            SP.PiCond(i,:) = p;
        end
        clear p

    else
        %-Full permutation test :
        % Build up exhaustive matrix of permutations
        %===============================================================
        %-Compute permutations for a single exchangability block
        %---------------------------------------------------------------
        %-Initialise XblkPiCond & remaining numbers
        XblkPiCond = [];
        lef = [1:stats.numsubj]';
        %-Loop through numbers left to add to permutations, accumulating PiCond
        for i = stats.numsubj:-1:1
            %-Expand XblkPiCond & lef
            tmp = round(exp(gammaln(stats.numsubj+1)-gammaln(i+1)));
            Exp = meshgrid(1:tmp,1:i); Exp = Exp(:)';
            if ~isempty(XblkPiCond), XblkPiCond = XblkPiCond(:,Exp); end
            lef = lef(:,Exp);
            %-Work out sampling for lef
            tmp1 = round(exp(gammaln(stats.numsubj+1)-gammaln(i+1)));
            tmp2 = round(exp(gammaln(stats.numsubj+1)-gammaln(i)));
            sam = 1+rem(0:i*tmp1-1,i) + ([1:tmp2]-1)*i;
            %-Add samplings from lef to XblkPiCond
            XblkPiCond   = [XblkPiCond; lef(sam)];
            %-Delete sampled items from lef & condition size
            lef(sam) = [];
            tmp = round(exp(gammaln(stats.numsubj+1)-gammaln((i-1)+1)));
            lef = reshape(lef,(i-1),tmp);
            %NB:gamma(stats.numsubj+1)/gamma((i-1)+1) == size(XblkPiCond,2);
        end
        clear lef Exp sam i
        %-Reorientate so permutations are in rows
        XblkPiCond = XblkPiCond';
        SP.PiCond=XblkPiCond;
    end

    %-Check, condition and randomise PiCond
    %-----------------------------------------------------------------------
    %-Check PiConds sum within Xblks to sum to 1
    if ~all(all(sum(SP.PiCond,2)== (stats.numsubj+1)*stats.numsubj/2 ))
        error('Invalid PiCond computed!')
    end
    %-Convert to full permutations from permutations within blocks
    stats.numperm = size(SP.PiCond,1);
    %-Randomise order of PiConds (except first) to allow interim analysis
    rand('seed',sum(100*clock))	%-Initialise random number generator
    SP.PiCond=[SP.PiCond(1,:);SP.PiCond(randperm(stats.numperm-1)+1,:)];
    %-Check first permutation is null permutation
    if ~all(SP.PiCond(1,:)==[1:stats.numsubj])
        error('SP.PiCond(1,:)~=[1:stats.numsubj]'); 
    end    
    
end     % switch snpm model
SP.dosmoothvar = ( stats.usespatnormed && (SP.STAT~='R') );
SP.G = stats.covar;

%==========================================================================
function subj=loadbeam(pa)

subj = load(pa);
ismsm = ( isfield(subj,'beam') && isfield(subj,'MNIsubjmap') );
if ismsm
    msm = subj.MNIsubjmap;
end
subj=nut_beam_legacy_compatibility(subj); 
if ismsm
    subj.MNIsubjmap = msm;
end
