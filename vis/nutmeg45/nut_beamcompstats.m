function nut_beamstats(stats)
% NUT_BEAMSTATS: calculates averages and statistics across subjects. 
%
% Usage:  nut_beamstats(stats)
%
% stats    Structure array created with nut_timef_stats GUI.


dimxyz=[91 109 91];
numvox=prod(dimxyz);  

ftxt=sprintf('_f%s',stats.frqtxt);      % freq info in new filename
if stats.frqstats, ftxt=[ftxt 'stat']; end
if strcmp(stats.fsel,'mean'), ftxt=[ftxt 'm']; end
if strcmp(ftxt,'_f'), ftxt=''; end

if stats.tsstats                                    % time info in filename
    ttxt=sprintf('_t%sstat',stats.timetxt);
elseif ~isempty(stats.timetxt)
    ttxt=sprintf('_t%s',stats.timetxt);
else, ttxt='';
end
clear MNIsubjmap

fprintf('\n\nNUT_BEAMCOMPSTATS: Compares and averages condition pairs.\n')
for cc=1:size(stats.comps,1)
    fprintf('\n\nComparing conditions %d and %d...\n',stats.comps(cc,:))
    tmpnr=sprintf('%dvs%d',stats.comps(cc,:));
    avgfile=['s_beamtf' tmpnr ftxt ttxt '_avg.mat'];
    subjcompfile=cell(1,stats.numsubj);
    clear tmpnr
    
    subj1=load(stats.files{stats.subjsel(1),stats.comps(cc,1)});        % load first subject
        
    if (isempty(stats.frqband) || stats.frqstats || stats.tsstats)
        stats.frqband=1:size(subj1.beam.bands,1);
    end
    if (isempty(stats.timepts) || stats.frqstats || stats.tsstats)
        stats.timepts=1:length(subj1.beam.timepts);
    end
    switch stats.fsel
        case 'one', newfreqbins=length(stats.frqband);
        case 'mean', newfreqbins=1;
    end
        
    if ~exist(avgfile,'file')       % skip if average has already been calculated in earlier session 
        switch stats.fsel
            case 'one'
                beam.bands = subj1.beam.bands(stats.frqband,:);
            case 'mean'
                beam.bands = [subj1.beam.bands(stats.frqband(1),1) subj1.beam.bands(stats.frqband(end),2)];
        end
        beam.timepts = subj1.beam.timepts(stats.timepts);
        beam.timewindow = subj1.beam.timewindow(stats.timepts,:);
        beam.srate = subj1.beam.srate;

        beam.coreg.mripath=[fileparts(which('nutmeg.m')) filesep 'templates' filesep 'wNormT1neuroax.img'];            %which('T1.img'); % should find SPM's T1 MNI template
        beam.coreg.brainrender_path=[fileparts(which('nutmeg.m')) filesep 'templates' filesep 'render_wNormT1neuroax_seg1.mat'];
        beam.coreg.meg2mri_tfm=eye(4);
        beam.coreg.orientation=1;
 
        MNItfm = [2 0 0 -92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];
        beam.voxels=nut_coordgrid(-90:2:90,-126:2:90,-72:2:108);
        MNIgrid=nut_coordtfm(beam.voxels,inv(MNItfm));

        beam.s{1} = zeros(numvox,length(stats.timepts),newfreqbins);
        MNIdeadvoxels = [];
        
        %------------------------------
        % Voxels common to both conditions and all subjects
        %------------------------------
        fprintf('\nFinding voxels common to both conditions all subjects, and calculating average...\n')
        for ii=1:stats.numsubj     % loop over datasets
            fprintf('\tProcessing subject %d...\n',stats.subjsel(ii))
            
            if(ii>1) % already loaded if ii = 1
                subj1=load(stats.files{stats.subjsel(ii),stats.comps(cc,1)});
            end
            subj2=load(stats.files{stats.subjsel(ii),stats.comps(cc,2)});
            
            subj1.beam.s{1}=subj1.beam.s{1}(:,stats.timepts,stats.frqband);
            subj2.beam.s{1}=subj2.beam.s{1}(:,stats.timepts,stats.frqband);
            
            if (~stats.frqstats && strcmp(stats.fsel,'mean'))
                subj1.beam.s{1}=mean(subj1.beam.s{1},3);
                subj2.beam.s{1}=mean(subj2.beam.s{1},3);
            end               
            
            % add NaN entries to end of subject beam, and remap missing MNI voxels to it   
            %subj1.beam.s{1}(end+1,:,:)=nan(1,length(stats.timepts),newfreqbins);
            %subj2.beam.s{1}(end+1,:,:)=nan(1,length(stats.timepts),newfreqbins);
            
            if isfield(subj1,'MNIsubjmap')
                MNIsubjmap=subj1.MNIsubjmap;
                subj1=rmfield(subj1,'MNIsubjmap');
                f=find(MNIsubjmap==size(subj1.beam.s{1},1)+1);      % this is a legacy of old algorithm
                if ~isempty(f), MNIsubjmap(f)=0; end                % set to 0
                clear f
            else
                voxelgrid = nut_coordtfm(subj1.beam.voxels,inv(MNItfm));
                [tmp,MNIsubjmap]=ismember(MNIgrid,voxelgrid,'rows');
                clear tmp voxelgrid
                %MNIsubjmap(find(MNIsubjmap==0))=size(subj1.beam.s{1},1);
                eval(['save ' stats.files{stats.subjsel(ii),stats.comps(cc,1)} ' MNIsubjmap -append']);
            end
            MNIsubjmap1=MNIsubjmap; clear MNIsubjmap
            deadvoxels = find(MNIsubjmap1==0);
            
            if isfield(subj2,'MNIsubjmap')
                MNIsubjmap=subj2.MNIsubjmap;
                subj2=rmfield(subj2,'MNIsubjmap');
                f=find(MNIsubjmap==size(subj2.beam.s{1},1)+1);      % this is a legacy of old algorithm
                if ~isempty(f), MNIsubjmap(f)=0; end                % set to 0
                clear f                
            else
                voxelgrid = nut_coordtfm(subj2.beam.voxels,inv(MNItfm));
                [tmp,MNIsubjmap]=ismember(MNIgrid,voxelgrid,'rows');
                clear tmp voxelgrid
                %MNIsubjmap(find(MNIsubjmap==0))=size(subj2.beam.s{1},1);
                eval(['save ' stats.files{stats.subjsel(ii),stats.comps(cc,2)} ' MNIsubjmap -append']);
            end
            MNIsubjmap2=MNIsubjmap; clear MNIsubjmap
            deadvoxels = union(deadvoxels,find(MNIsubjmap2==0));
            MNIdeadvoxels = union(MNIdeadvoxels,deadvoxels);
            clear deadvoxels
            
            %------------------
            % Calculate average
            %------------------        
            %for jj=1:length(beam.s)  %ii num of orientations
            subj.beam.s{1}=zeros(numvox,length(stats.timepts),newfreqbins);
            for vv=find(MNIsubjmap1 & MNIsubjmap2)  % vv num of voxels
                %if(MNIsubjmap(vv)~=0)
                subj.beam.s{1}(vv,:,:) = subj1.beam.s{1}(MNIsubjmap1(vv),:,:) - subj2.beam.s{1}(MNIsubjmap2(vv),:,:);
                %end
            end
            %end    
            clear subj1 subj2 MNIsubjmap1 MNIsubjmap2
            dum=pwd; cd(fileparts(stats.files{stats.subjsel(ii),stats.comps(cc,1)})); cd ..; cd ..;
            subjcompfile{ii}=fullfile(pwd,sprintf('s_beamtf_subj%dcond%dvs%d%s%s.mat',stats.subjsel(ii),stats.comps(cc,:),ftxt,ttxt));
            save(subjcompfile{ii},'subj');
            cd(dum); clear dum
            beam.s{1} = beam.s{1} + subj.beam.s{1};
            clear subj       
        end     % subject loop
        clear voxelgrid MNIgrid MNItfm
    
        % remove voxels not present in all subjects or 0 at all time points 
        Q = find(any(beam.s{1}(:,:),2));    % Q : voxels non 0 in 91*109*91 MNI space
        Q = setdiff(Q,MNIdeadvoxels)';      %     and defined in all subjects       
        clear MNIdeadvoxels
        beam.s{1}=beam.s{1}(Q,:,:);
        beam.voxels=beam.voxels(Q,:);
        beam.MNIgoodvoxels=Q;
        beam.subjcompfile=subjcompfile;

        beam.s{1} = beam.s{1}./stats.numsubj;  
        beam.voxelsize = [2 2 2];   

        fprintf('\tSaving average as %s...\n',avgfile)
        save(avgfile,'beam');
    elseif ~isempty(strmatch(stats.type,{'snpm','FDR'}))                   % load average calculated in former session for stats
        fprintf('\nUsing average file %s generated in earlier session.\n',avgfile)
        load(avgfile);
        Q=beam.MNIgoodvoxels;
        subjcompfile=beam.subjcompfile;
    end
    
    %-----------------
    % Stats part: SnPM
    %-----------------
    switch stats.type
    case 'snpm'
        fprintf('\nPerforming SnPM across subjects...\n')      
        
        VOX=beam.voxelsize;   
        clear beam 
        numgvox=length(Q);
        if isequal(stats.numperm,2^stats.numsubj), stats.hPerms=1;
        else, stats.hPerms=0;
        end
        
        % create permutation matrix
        fprintf('\tCreating permutation matrix...\n')
        if stats.numsubj<12 | stats.hPerms
            PiCond=[];
            for i=0:stats.numsubj-1
                PiCond=[ones(2^i,1),PiCond;-ones(2^i,1),PiCond];
            end
            clear i
            if stats.hPerms    %-Only do half the work, if possible
                PiCond=PiCond(PiCond(:,1)==1,:);
                numperm=stats.numperm/2;            % note: numperm is number of permutations actually calculated, stats.numperm the total number
            else
                numperm=stats.numperm;
                tmp=randperm(size(PiCond,1));
                PiCond=PiCond(tmp(1:numperm),:);
                clear tmp
            end
        else
            numperm=stats.numperm;
            d       = numperm-1;
            tmp     = pow2(0:stats.numsubj-1)*ones(stats.numsubj,1);  % Include correctly labeled iCond

            while (d>0)
              tmp = union(tmp,floor(rand(1,d)*2^stats.numsubj));
              tmp(tmp==2^stats.numsubj) = [];  % This will almost never happen
              d   = numperm-length(tmp);
            end

            % randomize tmp before it is used to get PiCond
            rand_tmp=randperm(length(tmp));
            tmp=tmp(rand_tmp);
            PiCond = 2*rem(floor(tmp(:)*pow2(-(stats.numsubj-1):0)),2)-1;
            clear tmp d rand_tmp
        end

        % initialize here only if possible, otherwise safe memory!
        nP    = ones(newfreqbins,numgvox); 
        %T    = zeros(newfreqbins,numgvox);
        %MaxT = zeros(stats.numperm,2,newfreqbins);       % permutations*[max min]*frqbands
                
        for ff=1:newfreqbins           % loop through frequency bands if requested.
            if newfreqbins>1, fprintf('\n\tPROCESSING FREQENCY BAND %d\n',stats.frqband(ff)), end
            
            fprintf('\tLoading all datasets ')
            %if stats.tsstats
            %    fprintf('and calculating mean of statistically significant time points...\n')
            %else
            fprintf('and calculating sum of all time points...\n')
            %end
            X = zeros(stats.numsubj,numgvox);
            for ii=1:stats.numsubj    % loop over datasets
                load(subjcompfile{ii});
                if newfreqbins>1
                    subj.beam.s{1}=subj.beam.s{1}(:,:,stats.frqband(ff));
                end
                X(ii,:)=sum(subj.beam.s{1}(Q,:),2);
                clear subj
            end  % subject loop
                           
            %---------
            % Variance
            %---------
            fprintf('\tCalculating variance of each voxel...\n')
            HC = PiCond(1,:)';
            
            BETA  = pinv(HC)*X;      % X: numsubj*num_non0voxels
            ResSS = sum((X - HC*BETA).^2);

            %-------------------
            % Variance smoothing
            %-------------------
            fprintf('\tSmoothing variance with a 20x20x20 mm gaussian kernel...\n')
            vFWHM     = [20 20 20];
            SmResSS   = zeros(dimxyz);
            SmMask    = zeros(dimxyz);
            TmpVol    = zeros(dimxyz);
            TmpVol(Q) = ones(size(Q));  % Q: non0voxels
            spm_smooth(TmpVol,SmMask,vFWHM./VOX);
            TmpVol(Q) = ResSS;
            spm_smooth(TmpVol,SmResSS,vFWHM./VOX);

            ResSS     = SmResSS(Q)./SmMask(Q);

            %--------------------------
            % Calculate t-values
            %--------------------------
            fprintf('\tCalculating T-values...\n')
            
            T(ff,:) = BETA ./ sqrt( (ResSS/stats.numsubj) / (stats.numsubj-1) );
            
            MaxT(1,:,ff) = [max(T(ff,:)) -min(T(ff,:))];
            %MaxT(1,:) = max(T(1,:));
            
            %if stats.hPerms
            %    nP(ff,:) = nP(ff,:) + (T(ff,:)<0);
            %end

            %-------------
            % Permutations
            %-------------  
            fprintf('\tPerforming permutations...\n')
            TP   = zeros(1,size(BETA,2));
            for perm = 2:numperm
                % Calculate variance
                HC = PiCond(perm,:)';
                BETA  = pinv(HC)*X;
                ResSS = sum((X - HC*BETA).^2);

                % Smooth Variance
                TmpVol(Q) = ResSS;
                spm_smooth(TmpVol,SmResSS,vFWHM./VOX);
                ResSS     = SmResSS(Q)./SmMask(Q);

                % Calculate t-values
                TP(1,:) = BETA ./ sqrt( (ResSS/stats.numsubj) / (stats.numsubj-1) );
                MaxT(perm,:,ff) = [ max(TP(1,:)) -min(TP(1,:)) ];
                %MaxT(perm,1) = max(TP(1,:));

                if stats.hPerms
                    nP(ff,:) = nP(ff,:) + (TP>=T(ff,:)) + (-TP>=T(ff,:));   % NB: Worry if T=TP=0
                    % if STAT=='T', then TP, T >=0, so (-TP>=T) will be empty.
                else
                    nP(ff,:) = nP(ff,:) + (TP>=T(ff,:));
                end

            end             
            clear X TP ResSS HC SmResSS SmMask BETA TmpVol perm
            
        end     % freqband loop
        
        %------------
        % p values
        %------------
        if newfreqbins>1, fprintf('\n'), end
        fprintf('\tCalculating p-values...\n')
        nP = nP/stats.numperm;

        nP_pos=nP;
        nP_neg = 1+1/stats.numperm-nP;
        clear nP
        
        if newfreqbins>1 
            if stats.corr4multfreq
                fprintf('\tCorrecting for testing multiple frequency bands.\n')
                MaxT = repmat(max(MaxT,[],3),[1 1 newfreqbins]);          % correct for multiple frequency bands
            else
                fprintf('\tDoes not correct for testing in multiple frequency bands.\n')
            end
        end
        
        if stats.hPerms
            MaxT = cat( 1, MaxT, flipdim(flipdim(MaxT,2),1) );
            %MaxT = [ MaxT; flipud(fliplr(MaxT)) ];
        end
        tol = 1e-4;	% Tolerance for comparing real numbers

        cP_pos = zeros(newfreqbins,numgvox);    
        for t=1:stats.numperm    %MaxT_pos'

            for ff=1:newfreqbins
            %-FEW-corrected p is proportion of randomisation greater or
            % equal to statistic.
            %-Use a > b -tol rather than a >= b to avoid comparing
            % two reals for equality.
                cP_pos(ff,:) = cP_pos(ff,:) + (MaxT(t,1,ff) > T(ff,:)-tol);      % correct for multiple voxels
            end
        end
        cP_pos = cP_pos / stats.numperm;  

        cP_neg = zeros(newfreqbins,numgvox);
        for t = 1:stats.numperm    %MaxT_neg'
            for ff=1:newfreqbins
                cP_neg(ff,:) = cP_neg(ff,:) + (MaxT(t,2,ff) > -T(ff,:)-tol);
            end
        end
        cP_neg = cP_neg / stats.numperm;  

        clear MaxT tol

        %------------------------------------
        % Saving stats in s_beam structure
        %------------------------------------
        fprintf('\tSaving results to %s...\n',avgfile)
        load(avgfile); 
        
        beam.snpm.T=T;
        %if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact'))
        beam.snpm.p_uncorr_pos=nP_pos;
        beam.snpm.p_corr_pos=cP_pos;
        %end
        %if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact'))
        beam.snpm.p_uncorr_neg=nP_neg;
        beam.snpm.p_corr_neg=cP_neg;
        %end   
        save(avgfile,'beam');           % resave average file with stats data, but no masking of nonsignificant voxels
        
        if (0)                          % saving of masked data not necessary, since nut_timef_viewer can perform masking
            tmpidx=strfind(avgfile,'avg');
            if stats.tsstats, tmpt='tstat_'; else, tmpt=''; end
            if ~isempty(strmatch(stats.whichtail,{'pos','both'},'exact'))
                for ff=1:newfreqbins
                    beam.s{1}(find(cP_pos(ff,:)>0.05),:,ff)=0;
                end
                snpmfile=[avgfile(1:tmpidx-1) tmpt 'snpm_pos.mat'];
                save(snpmfile,'beam');      % save average with masking
            end
            if ~isempty(strmatch(stats.whichtail,{'neg','both'},'exact'))
                for ff=1:newfreqbins
                    beam.s{1}(find(cP_neg(ff,:)>0.05),:,ff)=0;
                end
                snpmfile=[avgfile(1:tmpidx-1) tmpt 'snpm_neg.mat'];
                save(snpmfile,'beam');    
            end
        end
        clear beam T nP* cP* avgfile snpmfile tmpidx tmpt
                
    %-----------------------------------
    % Stats: t-tests with FDR correction    
    %-----------------------------------    
    case 'FDR'
        clear beam 
        numgvox=length(Q);
        fprintf('\nPerforming t-tests with FDR correction for multiple testing...\n')      
        
        tval   = zeros(newfreqbins,numgvox);    % intialize
        for ff=1:newfreqbins           % loop through frequency bands if requested.
            if newfreqbins>1, fprintf('\n\tPROCESSING FREQENCY BAND %d\n',stats.frqband(ff)), end
            
            fprintf('\tLoading all datasets ')
            %if stats.tsstats
            %    fprintf('and calculating mean of statistically significant time points...\n')
            %else
            fprintf('and calculating sum of all time points...\n')
            %end 
            X = zeros(stats.numsubj,numgvox);
            for ii=1:stats.numsubj    % loop over datasets
                load(subjcompfile{ii});
                if newfreqbins>1
                    subj.beam.s{1}=subj.beam.s{1}(:,:,stats.frqband(ff));
                end
                X(ii,:)=sum(subj.beam.s{1}(Q,:),2);
                clear subj
            end  % subject loop
           
            %--------------------------
            % Calculate t-values
            %--------------------------        
            fprintf('\tCalculating t-values...\n')
            BETA = mean(X,1);
            var  = (X - repmat(BETA,[stats.numsubj 1])).^2; 
            var  = sum(var); 
            stddev = var / (stats.numsubj - 1); clear var X
            stddev = sqrt(stddev);
            
            for vv=1:size(BETA,1)
                tval(ff,vv) = BETA(1,vv) ./ ( stddev(1,vv) / sqrt(stats.numsubj) );
            end
            clear stddev BETA
            
        end     % frqband loop

        %--------------------------
        % Calculate p-values
        %--------------------------    
        if newfreqbins>1, fprintf('\n'), end
        fprintf('\tCalculating p-values...\n')
        switch stats.whichtail
        case 'pos'    
            uP = tcdf(-tval,stats.numsubj-1); 
        case 'both'
            uP = 2*tcdf(-abs(tval),stats.numsubj-1);
        case 'neg'    
            uP = tcdf(tval,stats.numsubj-1); 
        end
        %clear tval 

        %------------------------------------------
        % correct for multiple testing using 5% FDR
        %------------------------------------------
        fprintf('\tCorrecting for multiple testing using a 5%% FDR...\n')  
        if (newfreqbins>1 && stats.corr4multfreq)
            numtests=numgvox*newfreqbins;
        else
            numtests=numgvox;
        end
        j=1:numtests;
        C_N=sum(ones(1,numtests)./j);
        
        if (newfreqbins==1 || stats.corr4multfreq)
            p_vec=reshape(uP,[1 numtests]);
            [A, ix]=sort(p_vec);    clear p_vec
            B=ix*stats.FDRvox/(C_N*numtests);
            C=A-B;
            f=find(C<0);
            if ~isempty(f)
                d=max(f);
                cutoff=A(d);
            else
                cutoff=0;
            end

            cP=min((A*numtests*C_N)./ix,1);
            cP(ix)=cP;    
            cP=reshape(cP,[newfreqbins numgvox]);
        else
            for ff=1:newfreqbins
                p_vec = uP(ff,:);
                [A, ix]=sort(p_vec);    clear p_vec
                B=ix*stats.FDRvox/(C_N*numtests);
                C=A-B;
                f=find(C<0);
                if ~isempty(f)
                    d=max(f);
                    cutoff=A(d);
                else
                    cutoff=0;
                end

                cPtmp=min((A*numtests*C_N)./ix,1);
                cPtmp(ix)=cPtmp;
                cP(ff,:)=cPtmp;
                clear cPtmp
            end
        end
        
        clear A B C C_N d f ix j numtests

        fprintf('\tSaving results to %s...\n',avgfile)
        load(avgfile);
        beam.ttestFDR.tail=stats.whichtail;
        beam.ttestFDR.T=tval;
        beam.ttestFDR.p_uncorr=uP;
        beam.ttestFDR.p_corr=cP;
        beam.ttestFDR.cutoff=cutoff;
        save(avgfile,'beam');
        
        if (0)
            for ff=1:newfreqbins
                beam.s{1}(find(uP(ff,:)>cutoff),:,ff)=0;
            end
            tmpidx=strfind(avgfile,'avg');
            ttestfile=[avgfile(1:tmpidx-1) 'ttest_' stats.whichtail '.mat'];
            save(ttestfile,'beam');
        end
        clear uP cutoff tval cP
        clear beam ttestfile tmpidx avgfile
        
    otherwise
        clear beam Q
        
    end     % switch stats.type
end     % condition loop

fprintf('Done.\n')
