 function nut_beamstats_notnormed(stats)
% NUT_BEAMSTATS_NOTNORMED: calculates averages and snpm statistics across subjects. 
%  All subjects must have the same voxel coordinates.
%
%    Usage:  nut_beamstats_notnormed(stats)
% stats   Structure array created with nut_timef_stats GUI.
%
% Parts of the code are adapted from the SnPM3 toolbox:  http://www.fil.ion.ucl.ac.uk/spm/snpm
%

% structures used:
% ---------------
% stats:        settings made in nut_timef_stats graphical user interface (GUI).
% SP:           snpm parameters passed to nut_snpm
% beam:         average beamformer data and stat results for display in nut_timef_viewer
% subj:         beamformer data of single subjects.

error('This function has been integrated into nut_beamstats and is obsolete.')

fprintf('\n\nNUT_BEAMSTATS_NOTNORMED: Averaging and statistics tool.\n')


% Get voxelsize and prepare initial parameters
%---------------------------------------------

temp=loadbeam(stats.files{stats.subjnr(1,1),stats.condnr(1,1)});  % load first subject to get voxelsize

% beam.voxelsize = [5 5 5];
voxelsize = temp.voxelsize;
voxels=temp.voxels;        
numvox = size(voxels,1);
numfiles=size(stats.subjnr,2);

if (isempty(stats.frqband) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.frqband=1:size(temp.bands,1);
end
if (isempty(stats.timepts) || stats.frqstats || stats.tsstats || stats.docomp)
    stats.timepts=1:length(temp.timepts);
end
clear temp

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
        fprintf('\nCalculating average...\n')
        
        % Prepare beam structure
        %-----------------------
        subj=loadbeam(stats.files{stats.subjnr(nn,1),stats.condnr(nn,1)});        % load first subject

        if length(subj.s)>2, subj=nut_calc_beamratio(subj); end
        
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
        beam.voxels = voxels;    
        
        if stats.doavg
            beam.s{1} = zeros(numvox,length(stats.timepts),newfreqbins);
            %allsubj   = zeros(numvox,stats.numsubj,length(stats.timepts));
        end
        beam.coreg=subj.coreg;
       

        for ii=1:numfiles     % loop over datasets
            fprintf('\tProcessing subject %d',stats.subjnr(nn,ii))
            if ~isempty(strmatch(stats.type,{'anova_snpm','ttest2_snpm'})), fprintf(', condition %d...\n',stats.condnr(ii))
            else, fprintf('...\n')
            end
            
            if ii>1     % already loaded otherwise
                subj=loadbeam(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)});
        
                if ~isequal(subj.voxelsize,voxelsize), error('NUT_BEAMSTATS_NOTNORMED: All subjects and conditions must have the same voxel size!'), end
                if length(subj.s)>2, subj=nut_calc_beamratio(subj); end
            end
            
            subj.s{1}=subj.s{1}(:,stats.timepts,stats.frqband);
            
            if (~stats.frqstats && strcmp(stats.fsel,'mean'))
                subj.s{1}=mean(subj.s{1},3);
            end                      
              
            %------------------
            % Calculate average
            %------------------     
            if stats.doavg
                for jj=1:length(beam.s)  %ii num of orientations
                    beam.s{jj} = beam.s{jj} + subj.s{jj};
                end
            end
            clear subj
            
        end     % subject loop
    
        % Save data calculated so far
        %----------------------------
        if stats.doavg
            beam.s{1} = beam.s{1}./stats.numsubj;  
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
    end
    
    %---------------------------------------
    % Stats part: SnPM 
    %---------------------------------------
    if stats.dostat  % if snpm selected
        Q = 1:numvox;
        
        if isfield(beam,'snpm')
            if ( (isfield(beam.snpm,'p_corr_pos') && ~isempty(strmatch(stats.whichtail,{'pos','both'}))) || (isfield(beam.snpm,'p_corr_neg') && ~isempty(strmatch(stats.whichtail,{'neg','both'}))) )
                answer=questdlg('WARNING: Results of SnPM tests with the current hypothesis are already present in your current s_beam*.mat file. If you continue, these results will be overwritten and lost. Do you wish to continue?', ...
                    'Genie Question','Yes, go ahead','Oh no, stop please','Yes, go ahead');
                if strcmp(answer,'Oh no, stop please'), error('NUT_BEAMSTATS: User abort.'), end
            end
        end
        beam = struct('voxels',beam.voxels,'voxelsize',beam.voxelsize,'coreg',beam.coreg);
        numgvox=length(Q);      

        SP = make_permmatrix(stats);
                
        % Loop through time windows and frequency bands
        %----------------------------------------------    
        nP   = ones(numgvox,newtimepts,newfreqbins);        % initialize
        T    = zeros(numgvox,newtimepts,newfreqbins);
        MaxT = zeros(SP.numperm,2,newfreqbins,newtimepts);       % permutations*[max min]*frqbands*timepts
        
        for tt=1:newtimepts            % loop through time windows if requested.
            if newtimepts>1, fprintf('\n\tPROCESSING TIME WINDOW %d\n',stats.timepts(tt)), end
            for ff=1:newfreqbins           % loop through frequency bands if requested.
                if newfreqbins>1, fprintf('\n\tPROCESSING FREQENCY BAND %d\n',stats.frqband(ff)), end

                fprintf('\tLoading all datasets...\n')
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
                X = zeros(numfiles,numgvox);       % numfiles = num_subjects*num_conditions in case of ANOVA
                for ii=1:numfiles    % loop over datasets
                    subj = loadbeam(stats.files{stats.subjnr(nn,ii),stats.condnr(nn,ii)});                    
                    if length(subj.s)>2, subj=nut_calc_beamratio(subj); end
    
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
                    subj=rmfield(subj,'beam');   

                    X(ii,:)=ts;
    
                    clear ts subj
                end  % subject loop
                
                % SnPM main part (permutations)
                %------------------------------
                [T(:,tt,ff),MaxT(:,:,ff,tt),nP(:,tt,ff)]=nut_snpm_nosmooth(X,Q,voxelsize,SP);

            end     % freqband loop
        end     % timepts loop
        if newfreqbins>1 || newtimepts>1, fprintf('\n'), end

        % p values
        %------------
        fprintf('\tCalculating p-values...\n')
        nP = nP/stats.numperm;

        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            nP_pos=nP;
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT=='T'
            nP_neg = 1+1/stats.numperm-nP;
        end
        clear nP
        
        if newfreqbins>1 
            if stats.corr4multfreq
                fprintf('\tCorrecting for testing multiple frequency bands.\n')
                MaxT = repmat(max(MaxT,[],3),[1 1 newfreqbins 1]);          % correct for multiple frequency bands
            else
                fprintf('\tDoes not correct for testing multiple frequency bands.\n')
            end
        end
        if newtimepts>1 
            if stats.corr4multtime
                fprintf('\tCorrecting for testing multiple time windows.\n')
                MaxT = repmat(max(MaxT,[],4),[1 1 1 newtimepts]);          % correct for multiple time windows
            else
                fprintf('\tDoes not correct for testing multiple time windows.\n')
            end
        end
        
        if SP.hPerms
            MaxT = cat( 1, MaxT, flipdim(flipdim(MaxT,2),1) );
            %MaxT = [ MaxT; flipud(fliplr(MaxT)) ];
        end
        tol = 1e-4;	% Tolerance for comparing real numbers

        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            cP_pos = zeros(numgvox,newtimepts,newfreqbins);    
            for pp=1:stats.numperm    %MaxT_pos'
                for tt=1:newtimepts
                    for ff=1:newfreqbins
                        %-FEW-corrected p is proportion of randomisation greater or
                        % equal to statistic.
                        %-Use a > b -tol rather than a >= b to avoid comparing
                        % two reals for equality.
                        cP_pos(:,tt,ff) = cP_pos(:,tt,ff) + (MaxT(pp,1,ff,tt) > T(:,tt,ff)-tol);      % correct for multiple voxels
                    end
                end
            end
            cP_pos = cP_pos / stats.numperm;  
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT=='T'
            cP_neg = zeros(numgvox,newtimepts,newfreqbins);
            for pp = 1:stats.numperm    %MaxT_neg'
                for tt=1:newtimepts
                    for ff=1:newfreqbins
                        cP_neg(:,tt,ff) = cP_neg(:,tt,ff) + (MaxT(pp,2,ff,tt) > -T(:,tt,ff)-tol);
                    end
                end
            end
            cP_neg = cP_neg / stats.numperm;  
        end
        clear MaxT tol


        % Save stats in beam structure
        %------------------------------------
        if stats.doavg
            beam=loadbeam(avgfile); 
        else
            beam=loadbeam(statfile);
        end
        
        beam.snpm.T=T;
        if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact')) || SP.STAT=='F'
            beam.snpm.p_uncorr_pos=nP_pos;
            beam.snpm.p_corr_pos=cP_pos;
        end
        if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact')) && SP.STAT=='T'
            beam.snpm.p_uncorr_neg=nP_neg;
            beam.snpm.p_corr_neg=cP_neg;
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
            save(avgfile,'beam');           % resave average file with stats data, no masking of nonsignificant voxels
        end

        fprintf('\tSaving statistics as %s...\n',statfile)
        beam.s{1}=beam.snpm.T;
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

    rand('seed',sum(100*clock));        % Initialise random number generator
    SP.sHCform  = ['spm_DesMtx(SP.CovInt(SP.PiCond(perm,:))-mean(SP.CovInt),',...
                    '''C'',''CovInt(centred)'')'];
    SP.CONT     = [1];        
    SP.STAT     = 'T';
    answer=inputdlg('Enter name of workspace variable with behavioural data: ');
    SP.CovInt   = evalin('base',answer{1});
        
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

