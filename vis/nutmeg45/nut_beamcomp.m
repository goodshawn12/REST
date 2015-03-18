function stats=nut_beamcomp(stats)
% NUT_BEAMCOMP: calculates differences between condition pairs. 
%
% Usage:  stats=nut_beamcomp(stats)
%
% stats    Structure array created with nut_timef_stats GUI.


% Create common voxel grid
%-------------------------
subj1=load(stats.files{stats.subjsel(1),stats.comps(1,1)});  % load first subject to get voxelsize
subj1=nut_beam_legacy_compatibility(subj1);
voxelsize =subj1.voxelsize;
xx = unique(subj1.voxels(:,1))';     % Make sure the template grid matches the data voxel grid
yy = unique(subj1.voxels(:,2))';
zz = unique(subj1.voxels(:,3))';
xx = [ fliplr(xx(1)-voxelsize(1):-voxelsize(1): -90) xx xx(end)+voxelsize(1):voxelsize(1): 90 ];
yy = [ fliplr(yy(1)-voxelsize(2):-voxelsize(2):-130) yy yy(end)+voxelsize(2):voxelsize(2): 90 ];
zz = [ fliplr(zz(1)-voxelsize(3):-voxelsize(3): -80) zz zz(end)+voxelsize(3):voxelsize(3):110 ];
voxels=nut_coordgrid(xx,yy,zz);
MNItfm = [voxelsize(1) 0 0 xx(1); 0 voxelsize(2) 0 yy(1); 0 0 voxelsize(3) zz(1); 0 0 0 1];
MNIgrid=nut_coordtfm(voxels,inv(MNItfm));

useroi = isfield(subj1,'rois');
if useroi
    roideffile = subj1.R.roidef; 
    numvox = length(subj1.R.roilabel);
else
    numvox = size(voxels,1);
end

% Prepare some variables
%-----------------------
if ~isempty(stats.frqtxt)
    ftxt=sprintf('_f%s',stats.frqtxt);      % freq info in new filename
else
    ftxt=''; 
end
if ~isempty(stats.timetxt)                  % time info in filename
    ttxt=sprintf('_t%s',stats.timetxt);
else, ttxt='';
end
if isempty(stats.frqband) 
    stats.frqband=1:size(subj1.bands,1);
end
if isempty(stats.timepts) 
    stats.timepts=1:length(subj1.timepts);
end

subjcompfile=cell(max(stats.subjsel),size(stats.comps,1));
stats.MNIdeadvoxels=cell(1,size(stats.comps,1));
doask=true;

% Loop through condition pairs and subjects
%------------------------------------------
fprintf('\n\nNUT_BEAMCOMP: Calculates differences between condition pairs.\n')
for cc=1:size(stats.comps,1)
    fprintf('\n\nComparing conditions %d and %d...\n',stats.comps(cc,:))
    %tmpnr=sprintf('%dvs%d',stats.comps(cc,:));
    %avgfile=['s_beamtf' tmpnr ftxt ttxt '_avg.mat'];
    %clear tmpnr
%     if cc>1     % already loaded if first condition pair
%         subj1=load(stats.files{stats.subjsel(1),stats.comps(cc,1)});        % load first subject of condition pair
%         subj1=nut_beam_legacy_compatibility(subj1);
%         if ~isequal(subj1.voxelsize,voxelsize), error('NUT_BEAMCOMP: All subjects and conditions must have the same voxelsize!'), end
%     end
%     if useroi, subj1.s=subj1.rois; end
%     if length(subj1.s)>2,subj1=nut_calc_beamratio(subj1); end
        
    MNIdeadvoxels = [];

    % Voxels common to both conditions
    %------------------------------
    fprintf('\nFinding voxels common to both conditions, and calculating difference...\n')
    for ii=stats.subjsel     % loop over datasets
        fprintf('\tProcessing subject %d...\n',ii)
        
        % Determine where to put the contrast files (in parent directory of both conditions)
        p1 = fileparts(stats.files{ii,stats.comps(cc,1)}); p2 = fileparts(stats.files{ii,stats.comps(cc,2)});
        lp = min(length(p1),length(p2));
        commonpathidx = find( p1(1:lp) == p2(1:lp) ); 
        if ~isempty(commonpathidx)
            if ( isequal(length(commonpathidx),length(p1)) || strcmp(p1(commonpathidx(end)),filesep) ), pa=p1(commonpathidx);
            else pa=fileparts(p1(commonpathidx));
            end 
        else
            pa=p1;
        end   
        subjcompfile{ii,cc}=fullfile(pa,sprintf('s_beamtf%dvs%d_subj%d%s%s.mat',stats.comps(cc,:),ii,ftxt,ttxt));
        clear p1 p2 pa lp commonpathidx
        
        if ii==stats.subjsel(1), doskip = exist(subjcompfile{ii,cc},'file'); else doskip = doskip && exist(subjcompfile{ii,cc},'file'); end
        if doskip                                % ask if already done
            if doask
                answer=menu(sprintf('A file named %s already exists from a previous session. What would you like to do?',subjcompfile{ii,cc}), ...
                       'Recalculate the data and replace the existing file', ...
                       'Use the existing file', ...
                       'Always recalculate the data and replace the existing file', ...
                       'Always use existing file');
                pause(0.1);
                if answer<1, error('User abort.'), end
                doskip = (rem(answer,2)==0);
                doask  = (answer<3);
            end
        end
        if ~doskip
            if cc>1 || ii>1     % already loaded otherwise
                subj1=load(stats.files{ii,stats.comps(cc,1)});
                subj1=nut_beam_legacy_compatibility(subj1);                
                if ~isequal(subj1.voxelsize,voxelsize), error('NUT_BEAMCOMP: All subjects and conditions must have the same voxelsize!'), end
            end
            if useroi, subj1.s=subj1.rois; end
            if length(subj1.s)>1, subj1=nut_calc_beamratio(subj1); end
            subj1.s{1}=subj1.s{1}(:,stats.timepts,stats.frqband);
            
            voxelgrid = nut_coordtfm(subj1.voxels,inv(MNItfm));
            [tmp,MNIsubjmap1]=ismember(MNIgrid,voxelgrid,'rows');
            clear tmp voxelgrid
            %MNIsubjmap(find(MNIsubjmap==0))=size(subj1.s{1},1);
            %eval(['save ''' stats.files{ii,stats.comps(cc,1)} ''' MNIsubjmap -append']);
            %end
            subj1=rmfield(subj1,'voxels');
            %MNIsubjmap1=MNIsubjmap; clear MNIsubjmap
            deadvoxels = find(MNIsubjmap1==0);

            subj2=load(stats.files{ii,stats.comps(cc,2)});
            subj2=nut_beam_legacy_compatibility(subj2);
            if ~isequal(subj2.voxelsize,voxelsize), error('NUT_BEAMCOMP: All subjects and conditions must have the same voxel size!'), end
            if useroi, subj2.s=subj2.rois; end
            if length(subj2.s)>1, subj2=nut_calc_beamratio(subj2); end
            subj2.s{1}=subj2.s{1}(:,stats.timepts,stats.frqband);
            
            voxelgrid = nut_coordtfm(subj2.voxels,inv(MNItfm));
            [tmp,MNIsubjmap2]=ismember(MNIgrid,voxelgrid,'rows');
            clear tmp voxelgrid
            %MNIsubjmap(find(MNIsubjmap==0))=size(subj2.s{1},1);
            %eval(['save ''' stats.files{ii,stats.comps(cc,2)} ''' MNIsubjmap -append']);
            %end
            subj2=rmfield(subj2,{'voxels','coreg'});
            %MNIsubjmap2=MNIsubjmap; clear MNIsubjmap
            deadvoxels = union(deadvoxels,find(MNIsubjmap2==0));
            save([subjcompfile{ii,cc}(1:end-4) '_dv.mat'],'deadvoxels');      
            MNIdeadvoxels = union(MNIdeadvoxels,deadvoxels);
            clear deadvoxels 

            % Calculate Difference
            %------------------          
            if useroi
                [beam.R.goodroi,i1,i2] = intersect(subj1.R.goodroi,subj2.R.goodroi);
                beam.R.roidef = subj1.R.roidef;
                beam.R.roilabel = subj1.R.roilabel;
                beam.rois = { subj1.s{1}(i1,:,:) - subj2.s{1}(i2,:,:) };
                clear i1 i2
            else
                % save as much memory as possible, in expense of speed!
                vv=find(MNIsubjmap1 & MNIsubjmap2);
                MNIsubjmap1=MNIsubjmap1(vv);
                MNIsubjmap2=MNIsubjmap2(vv);

                B = zeros(length(vv),length(stats.timepts),length(stats.frqband));
                for kk=1:length(vv)
                    %if(MNIsubjmap(vv)~=0)
                    %subj1.s{1}(MNIsubjmap1(kk),:,:) = subj1.s{1}(MNIsubjmap1(kk),:,:) - subj2.s{1}(MNIsubjmap2(kk),:,:);
                    B(kk,:,:) = subj1.s{1}(MNIsubjmap1(kk),:,:) - subj2.s{1}(MNIsubjmap2(kk),:,:);
                    %end
                end  
                %subj1.s{1}=subj1.s{1}(unique(MNIsubjmap1),:,:);
                clear subj2 MNIsubjmap1 MNIsubjmap2 kk

                % bring contrast map to MNI space
                %----------------------------------
                beam.s{1}=zeros(numvox,length(stats.timepts),length(stats.frqband));
                beam.s{1}(vv,:,:)=B; clear B vv
            end
            beam.voxelsize = voxelsize;
            beam.voxels = voxels;
            beam.bands = subj1.bands(stats.frqband,:);
            beam.timepts = subj1.timepts(stats.timepts);
            beam.timewindow = subj1.timewindow(stats.timepts,:);
            beam.srate = subj1.srate;
            beam.coreg = subj1.coreg;

            clear subj1
            
            save(subjcompfile{ii,cc},'beam');
            clear beam       
        else
            clear subj1
            load([subjcompfile{ii,cc}(1:end-4) '_dv.mat']);
            MNIdeadvoxels = union(MNIdeadvoxels,deadvoxels);
            clear beam
        end
    end     % subject loop
    stats.MNIdeadvoxels{cc}=MNIdeadvoxels;          % store voxels not present in all subjects and conditions
end     % comps loop

stats.files=subjcompfile;                       % Update settings for subsequent analyses
stats.condsel=1:stats.numtests;
stats.subjnr=repmat(stats.subjsel,[stats.numtests 1]);
stats.condnr=repmat(stats.condsel',[1 stats.numsubj]);

