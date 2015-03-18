function stats=nut_ttest_freqband(stats)
% NUT_TTEST_FREQBAND: performs t-tests with the time-series of each frequency band
% and each voxel. For each voxel, the frequency bands with significant changes over 
% time are then selected and averaged. 
% Usage:  stats=nut_ttest_freqband(stats)
%         stats    Structure array created with nut_timef_stats GUI.

fprintf('\n\nNUT_TTEST_FREQBAND: Finds significant frequency bands and time points\n')
    fprintf('                    at each voxel and for each subject.\n')

if stats.frqstats
    tmpf=sprintf('_f%sstat',stats.frqtxt);     % freq info in new filename
elseif ~isempty(stats.frqtxt)
    tmpf=sprintf('_f%s',stats.frqtxt);     
else
    tmpf='';
end
if (strcmp(stats.fsel,'mean') && ~isempty(tmpf))
    tmpf=[tmpf 'm'];
end
if ~isempty(stats.timetxt)
    tmpt=sprintf('_t%s',stats.timetxt);     % time info in new filename
else
    tmpt='';
end

clear global beam       % to avoid problems when loading local beam structure

newfil=cell(size(stats.files));
for c=stats.condsel
    if length(stats.condsel)>1, fprintf('\n\nAnalyzing condition %d...\n',c), end
    for s=stats.subjsel
        if ~isempty(find(stats.subjnr==s & stats.condnr==c))
            fprintf('\nAnalyzing subject %d...\n',s)

            if stats.frqstats
                [fpath,fname,fext]=fileparts(stats.files{s,c});
                if ~isempty(stats.frqtxt)
                    f=strfind(fname,sprintf('_f%s',stats.frqtxt));      % do not put freq info twice in filename
                    if ~isempty(f), fname=fname(1:f-1); end
                end
                if ~isempty(stats.timetxt)
                    f=strfind(fname,sprintf('_t%s',stats.timetxt));     % do not put time info twice in filename
                    if ~isempty(f), fname=fname(1:f-1); end
                end
                newfil{s,c}=fullfile(fpath,[fname tmpf tmpt '_' stats.whichtail '.mat']);
                clear f*
                isnewfile=1;
            elseif (~isempty(tmpf) || ~isempty(tmpt))
                [fpath,fname,fext]=fileparts(stats.files{s,c});
                if ~isempty(stats.frqtxt)
                    f=strfind(fname,sprintf('_f%s',stats.frqtxt));
                    if ~isempty(f), fname=fname(1:f-1); end
                end
                if ~isempty(stats.timetxt)
                    f=strfind(fname,sprintf('_t%s',stats.timetxt));
                    if ~isempty(f), fname=fname(1:f-1); end
                end
                newfil{s,c}=fullfile(fpath,[fname tmpf tmpt '.mat']);
                clear f*    
                isnewfile=1;
            else
                isnewfile=0;
            end

            if (isnewfile && ~stats.docomp && exist(newfil{s,c},'file'))        % Skip if file was already created in earlier session
                [fpath,ffile]=fileparts(newfil{s,c});                           % Do not skip if stats.docomp, otherwise calculation is wrong
                fprintf('\tUsing file %s created in earlier session.\n',ffile)
                clear fpath ffile
            else

            fprintf('\tLoading dataset...\n')
            beam=load(stats.files{s,c});
            beam=nut_beam_legacy_compatibility(beam);
            if length(beam.s)>2, beam=nut_calc_beamratio(beam); end

            if stats.docomp
                beam.s{1}(stats.MNIdeadvoxels{c},:,:)=[];
                beam.voxels(stats.MNIdeadvoxels{c},:)=[];
                numvox=size(beam.s{1},1);
            else
                [numvox,numtim,numfrq]=size(beam.s{1});
                if isempty(stats.frqband)
                    stats.frqband=1:numfrq;
                end
                if isempty(stats.timepts)
                    stats.timepts=1:numtim;
                end
                beam.s{1}=beam.s{1}(:,stats.timepts,stats.frqband);
                beam.bands=beam.bands(stats.frqband,:);
                beam.timepts=beam.timepts(stats.timepts);
                beam.timewindow=beam.timewindow(stats.timepts,:);
            end
            numfrq=length(stats.frqband);
            numtim=length(stats.timepts);

            % test significance of each time point at each voxel
            if stats.tsstats
                warning('NUT_TTEST_FREQBAND: Warning, this option is not fully tested, and generally not recommended!')
                fprintf('\tTesting significance of TIME POINTS with t-tests...\n')
                if ~isfield(beam,'goodtimepts')
                    if numfrq<4
                        error('Dataset does not have enough frequency bands to perform t-tests.')
                    end

                    tval=zeros(numvox,numtim);
                    for tt=1:numtim
                        tval(:,tt) = mean(beam.s{1}(:,tt,:),3) ./ ( std(beam.s{1}(:,tt,:),[],3)./sqrt(numfrq) + 0.0000001 );
                    end

                    switch stats.whichtail
                    case 'pos'  
                        p = ones(size(tval));
                        for tt=1:size(tval,2)
                            p(:,tt) = tcdf(-tval(:,tt),numfrq-1); 
                        end
                    case 'neg'
                        p = ones(size(tval));
                        for tt=1:size(tval,2)
                            p(:,tt) = tcdf(tval(:,tt),numfrq-1);
                        end
                    otherwise
                        error('Invalid tail for nut_ttest_freqband.')
                    end
                    clear tval

                    beam.goodtimepts=cell(numvox,1);
                    if numtim>1         % correct for multiple testing over several frequency bands using 5% FDR
                        fprintf('\tFinding significant time points...\n')
                        cutoff=0.05;
                        if stats.corr4multtime
                            fprintf('\tApplying Bonferroni correction for testing multiple time windows.\n')
                            cutoff=cutoff/numtim;
                        end

                        for k=1:numvox
                            beam.goodtimepts{k}=find(p(k,:)<cutoff);   
                        end
                        clear cutoff

                    else
                        fprintf('\tFinding voxels with significant frequency bands at tested time point...\n')
                        beam.goodtimepts(find(p<0.05))={1};
                    end

                    if ~stats.frqstats      % Save in beam structure of original file if no freq stats performed (otherwise save later)
                        fprintf('\tSaving...\n')
                        if isnewfile
                            save(newfil{s,c},'beam')
                        else
                            save(stats.files{s,c},'beam')
                        end
                    end
                else
                    fprintf('\tUsing values calculated in earlier session.\n',ffile)
                end
            end     % if stats.tsstats

            % test significance of each frequency band at each voxel
            if stats.frqstats
                fprintf('\tTesting significance of FREQUENCY BANDS with t-tests...\n')

                if numtim<4
                    error('Dataset does not have enough time points to perform t-tests.')
                end

                tval=zeros(numvox,numfrq);
                for ff=1:numfrq
                    tval(:,ff) = mean(beam.s{1}(:,:,ff),2) ./ ( std(beam.s{1}(:,:,ff),[],2)./sqrt(numtim) + 0.0000001 );
                end

                switch stats.whichtail
                case 'pos'    
                    p = tcdf(-tval,numtim-1); 
                case 'neg'
                    p = tcdf(tval,numtim-1);
                otherwise
                    error('Invalid tail for nut_ttest_freqband.')
                end
                clear tval

                if numfrq>1         % correct for multiple testing over several frequency bands using FDR
                   switch stats.fsel
                        case 'one'
                            TS=zeros(numvox,numtim,numfrq);    % initialize new time series
                            fprintf('\tFinding significant frequency bands...\n')
                        case 'mean'
                            TS=zeros(numvox,numtim);
                            fprintf('\tAveraging significant frequency bands...\n')
                   end
                   cutoff=0.05;
                   if stats.corr4multfreq
                       fprintf('\tApplying Bonferroni correction for testing multiple frequency bands.\n')
                       cutoff=cutoff/numfrq;
                   end
                   for k=1:numvox
                       f=find(p(k,:)<cutoff);
                       if ~isempty(f)
                           switch stats.fsel
                               case 'one'
                                   TS(k,:,f)=beam.s{1}(k,:,f);
                               case 'mean'
                                   TS(k,:)=mean(beam.s{1}(k,:,f),3);
                           end
                       end
                   end

                   beam.s{1}=TS; clear TS cutoff
                else
                    fprintf('\tFinding voxels with significant time series at tested frequency band...\n')
                    beam.s{1}(find(p>=0.05),:)=0;
                end

                %stats.timepts=1:length(beam.timepts);
                if strcmp(stats.fsel,'mean')
                    beam.bands=[beam.bands(1,1) beam.bands(end,2)];
                end

                fprintf('\tSaving data...\n')
                save(newfil{s,c},'-struct','beam')        % always isnewfile
                    %clear beam 
            end     % if stats.frqstats

            end     % if skip
        end  % if ~isempty
    end     % for subj
end     % for condition

if (nargout>0 && isnewfile)
    stats.files=newfil;
end
