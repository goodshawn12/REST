function nut_render_meg(rendfile)
% NUT_RENDER_MEG(RENDFILE)
%
% Optional to put in RENDFILE (see spm_render).  If not specified, will be
%         prompted for.  The RENDFILE should be created first using spm.
%
% Uses global ST
%
% only makes sense with activations viewed in neurological format
% (otherwise will be flipped)

global st rivets
if(~isfield(st.vols{1},'blobs'))
    [dat(1:3).mat]=deal(eye(4));
    [dat(1:3).t]=deal(zeros(4,1));
    [dat(1:3).XYZ]=deal(ones(3,4));
    [dat(1:3).dim]=deal([4 1 1]);
    if nargin < 1,
        rendfile = spm_get(1,'render*.mat','Render file',fullfile(spm('Dir'),'rend'));
    end;
    spm_render(dat,1,rendfile);
    return
    
else 
    [dat(1:3).mat]=deal(st.vols{1}.blobs{1}.mat);
%     dat(2).mat=st.vols{1}.blobs{1}.mat;
%     dat(3).mat=st.vols{1}.blobs{1}.mat;
    xdim=size(st.vols{1}.blobs{1}.vol,1);
    ydim=size(st.vols{1}.blobs{1}.vol,2);
    zdim=size(st.vols{1}.blobs{1}.vol,3);
    numblobvoxels=xdim*ydim*zdim;
    [dat(1:3).dim]=deal([xdim,ydim,zdim]);
%     dat(2).dim=[xdim,ydim,zdim];
%     dat(3).dim=[xdim,ydim,zdim];

    %[i,j,k]=ind2sub(size(st.vols{1}.blobs{1}.vol),find(st.vols{1}.blobs{1}.vol>max(rivets.threshold(:,1))));
    [i,j,k]=ind2sub(size(st.vols{1}.blobs{1}.vol),find(st.vols{1}.blobs{1}.vol>0));     % Thresholding has already been done before and need not be repeated here
    dat(1).XYZ = [i';j';k'];
    numpos=length(i);
    % dat(2).XYZ = dat(1).XYZ;
    %[i,j,k]=ind2sub(size(st.vols{1}.blobs{1}.vol),find(st.vols{1}.blobs{1}.vol<min(rivets.threshold(:,2))));
    [i,j,k]=ind2sub(size(st.vols{1}.blobs{1}.vol),find(st.vols{1}.blobs{1}.vol<0));     % Thresholding has already been done before and need not be repeated here
    dat(3).XYZ = [i';j';k'];
    numneg=length(i);
    clear i j k;
end
% a=[];
% b=[];
% ind=1;
% for kk=1:zdim,
%     [a,b]=find(st.vols{1}.blobs{1}.vol(:,:,kk)>0);  %anything above 0 will be the threshold specified in current view
%     if ~isempty(a)
%         for ii=1:length(a),
%             dat.XYZ(:,ind)=[a(ii) b(ii) kk]';
%             ind=ind+1;
%         end
%     end
%     clear a b
% end

if numpos>0
    ind=sub2ind(size(st.vols{1}.blobs{1}.vol),dat(1).XYZ(1,:),dat(1).XYZ(2,:),dat(1).XYZ(3,:));
    dat(1).t=st.vols{1}.blobs{1}.vol(ind)'; % was *1.5
end
if numpos<4
% if spm_render has less than 4 voxels, it goes psycho. Therefore, we fill up with 0s
% here
    dat(1).t(numpos+1:4,1) = 0;
    dat(1).XYZ(1:3,numpos+1:4)=1;
end

dat(2).t=zeros(size(dat(1).t));
dat(2).XYZ=dat(1).XYZ;

if numneg>0
    ind=sub2ind(size(st.vols{1}.blobs{1}.vol),dat(3).XYZ(1,:),dat(3).XYZ(2,:),dat(3).XYZ(3,:));
    dat(3).t=abs(st.vols{1}.blobs{1}.vol(ind))';    
end
if numneg<4
% if spm_render has less than 4 voxels, it goes psycho. Therefore, we fill up with 0s
% here
    dat(3).t(numneg+1:4,1) = 0;
    dat(3).XYZ(1:3,numneg+1:4)=1;
end
% 	for ii=1:size(dat.XYZ,2)
%         dat.t(ii)=st.vols{1}.blobs{1}.vol(dat.XYZ(1,ii),dat.XYZ(2,ii),dat.XYZ(3,ii));
% 	end
% 	dat.t=dat.t';

if nargin < 1
    if strcmp(spm('ver'),'SPM2')
        rendfile = spm_get(1,'render*.mat','Render file',fullfile(spm('Dir'),'rend'));
    else
        rendfile = [];
    end
else
    if strncmp(spm('ver'),'SPM8',4)     % SPM8 or 8b
        global prevrend
        if( ~isstruct(prevrend) || isempty(prevrend.rendfile) )
            prevrend = struct('rendfile',rendfile,'brt',[],'col',[]);
        end
    end
end;

spm_render(dat,1,rendfile);



