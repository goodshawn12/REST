function VOIvoxels = nut_select_VOI(asksave)
% NUT_SELECT_VOI
%
% This function is called when "Select VOI" is pressed.
%
% Currently requires image processing toolbox.
% TODO: rbbox can be used to draw rectangular ROIs (esp. for people who don't
% have this toolbox.)
%
% Uses globals st and nuts

% note that spm_orthviews is not singles friendly. it will crash and take
% matlab with it, like an Indian mom with an unmarried 30-year-old kid.
                                                                                
global nuts st ndefaults

if nargin<1, asksave=true; end

if(~exist('roipoly.m','file')) % check for image processing toolbox
    warning('Image Processing Toolbox not found -- graphical selection unavailable');
    ndefaults.voi.selectmethod = 3;
end

switch(ndefaults.voi.selectmethod)
    case 1
        use_roipoly = true;
    case 2  % use drawrect instead
        use_roipoly = false;
    case 3
        prompt={'[xmin xmax]'
                '[ymin ymax]'
                '[zmin zmax]'};
        name='Choose VOI bounds (MRI coords)...';
        if(~exist('answer','var'))
            answer=inputdlg(prompt,name); %,numlines,defaultanswer);
        end

        if(isempty(answer))
            return; % no input, let's blow this taco stand
        end
        
        for i=1:3
            answer{i}=str2num(answer{i});
        end
        [xlims,ylims,zlims] = deal(answer{:});
        VOIvoxels = nut_coordgen([xlims ylims zlims]);

        %%%%%%% display chosen VOI
        blobcenter = mean(VOIvoxels);

        % spm blobs need to have positive coordinate values
        translation_tfm = [ 1 0 0 -abs(st.bb(1,1))
                            0 1 0 -abs(st.bb(1,2))
                            0 0 1 -abs(st.bb(1,3))
                            0 0 0                1 ];
        blobvoxels = nut_coordtfm(VOIvoxels,inv(translation_tfm));

        spm_image('init',st.vols{1}.fname); % reload structural MRI; clears VOI emphasis from nut_view_scan_region.m

        spm_orthviews('rmblobs',1)
        spm_orthviews('Reposition',blobcenter);

        %  so that the coloured blobs are not opaque, a contrast within the blob
		%  needs to exist; hence arbitrarily setting the last voxel to value=4
		%  while the rest remain value=1
        spm_orthviews('addcolouredblobs',1,blobvoxels',ones(size(blobvoxels,1),1),translation_tfm,[.3 1 .2]);
        st.vols{1}.blobs{1}.max=4;  % force desired color scale
        st.vols{1}.blobs{1}.min=0;
	    spm_orthviews('redraw');
        
        % VOIvoxels saved in MEG coords
        VOIvoxels = nut_mri2meg(VOIvoxels);

        [filename, pathname] = uiputfile('*.mat', 'Save VOI...(OPTIONAL)','voxels.mat');  %optional as long as saving VOIvoxels to nuts
        if ~(isequal(filename,0) | isequal(pathname,0))
            save(fullfile(pathname, filename),'VOIvoxels');
        end

        % remove any lead field since we've obviously made it invalid by loading a new VOI
        if(isfield(nuts,'Lp'))
            nuts = rmfield(nuts,'Lp');
        end

        nut_enabler;
        return
    otherwise
        errordlg('Unsupported VOI selection method specified in ndefaults.');       
end



spm_orthviews('rmblobs',1);

 %%%% NOTE -- SPM axes are in mm, with indices beginning at lower left of image

 %%% TODO: doing something like this can drastically reduce the size of
 %%% VOIvoxels..
 % [a,x,y] = roipoly;
 % aa=roipoly(get(st.vols{1}.ax{1}.d,'Cdata'),x,y)
 % this results in aa = a...
 
% get VOI dims from st.bb (bounding box)
VOIdim = round(diff(st.bb)+1);

translation_tfm = [ 1 0 0 -abs(st.bb(1,1))
                    0 1 0 -abs(st.bb(1,2))
                    0 0 1 -abs(st.bb(1,3))
                    0 0 0                1 ];
 
axes_properties = {'LineWidth','XColor','YColor'};
outlined_axes_values = {3, [1 0 0], [1 0 0]};

% axial view (x,y)
axes(st.vols{1}.ax{1}.ax);
restore_axes_values = get(st.vols{1}.ax{1}.ax,axes_properties);
set(st.vols{1}.ax{1}.ax, axes_properties,outlined_axes_values);
texthandle=uicontrol('Style','text','BackgroundColor',[1 1 1],'String','Select axial cross-section','Units','normalized','Position',[.6 .55 .3 .1],'FontWeight','bold','FontSize',16);

if(use_roipoly)
    ROI_ax = single(roipoly');   %ends up in mm coordinates!

    VOI=repmat(ROI_ax,[1 1 VOIdim(3)]);
    clear ROI_ax
    [VOIx,VOIy,VOIz]=ind2sub(size(VOI),single(find(VOI)));
else
    voirect = single(getrect);
    rectx = round(voirect(1):(voirect(1)+voirect(3)));
    recty = round(voirect(2):(voirect(2)+voirect(4)));
    [VOIx,VOIy,VOIz] = meshgrid(rectx, recty, 1:VOIdim(3));
end
VOIvoxels = [VOIx(:) VOIy(:) VOIz(:)];clear VOIx VOIy VOIz;

set(st.vols{1}.ax{1}.ax,axes_properties,restore_axes_values);
%  so that the coloured blobs are not opaque, a contrast within the blob
%  needs to exist; hence we arbitrarily set the max value to 4
cursorpos = spm_orthviews('Pos');
spm_orthviews('rmblobs',1)
if ndefaults.voi.reposition
    blobcenter = double(nut_coordtfm(mean(VOIvoxels),translation_tfm));
    spm_orthviews('Reposition',[blobcenter(1) blobcenter(2) cursorpos(3)]);
end
spm_orthviews('addcolouredblobs',1,VOIvoxels',ones(size(VOIvoxels,1),1),translation_tfm,[.3 1 .2]);
st.vols{1}.blobs{1}.max=4;  % force desired color scale
st.vols{1}.blobs{1}.min=0;
spm_orthviews('redraw');
clear VOIvoxels


% coronal view (x,z)
axes(st.vols{1}.ax{2}.ax);
restore_axes_values = get(st.vols{1}.ax{2}.ax,axes_properties);
set(st.vols{1}.ax{2}.ax, axes_properties,outlined_axes_values);
set(texthandle,'String','Select coronal cross-section');
if(use_roipoly)
    ROI_co = single(roipoly');
    VOI_co=reshape(ROI_co,[size(ROI_co,1) 1 size(ROI_co,2)]);
    VOI=VOI.*repmat(VOI_co,[1 VOIdim(2) 1]);
    clear ROI_co VOI_co
    [VOIx,VOIy,VOIz]=ind2sub(size(VOI),single(find(VOI)));
else
    voirect = single(getrect);
    rectx = intersect(rectx,round(voirect(1):(voirect(1)+voirect(3))));
    rectz = round(voirect(2):(voirect(2)+voirect(4)));

    [VOIx,VOIy,VOIz] = meshgrid(rectx, recty, rectz);
end
VOIvoxels = [VOIx(:) VOIy(:) VOIz(:)]; clear VOIx VOIy VOIz;

set(st.vols{1}.ax{2}.ax,axes_properties,restore_axes_values);
cursorpos = spm_orthviews('Pos');
spm_orthviews('rmblobs',1)
if ndefaults.voi.reposition
    blobcenter = double(nut_coordtfm(mean(VOIvoxels),translation_tfm));
    spm_orthviews('Reposition',[blobcenter(1) cursorpos(2) blobcenter(3)]);
end
spm_orthviews('addcolouredblobs',1,VOIvoxels',ones(size(VOIvoxels,1),1),translation_tfm,[.3 1 .2]);
st.vols{1}.blobs{1}.max=4;  % force desired color scale
st.vols{1}.blobs{1}.min=0;
spm_orthviews('redraw');
clear VOIvoxels


% sagittal view (y,z)
axes(st.vols{1}.ax{3}.ax);
restore_axes_values = get(st.vols{1}.ax{3}.ax,axes_properties);
set(st.vols{1}.ax{3}.ax, axes_properties,outlined_axes_values);
set(texthandle,'String','Select sagittal cross-section');

if(use_roipoly)
    ROI_sa = single(roipoly');
    ROI_sa = flipud(ROI_sa); % y is flipped in this view
    VOI_sa=reshape(ROI_sa,[1 size(ROI_sa,1) size(ROI_sa,2)]);
    %pack; % just to be nice to microsoft's victims. poor fools don't have decent memory management.
    clear ROI_sa
    VOI=VOI.*repmat(VOI_sa,[VOIdim(1) 1 1]);
    clear VOI_sa
    [VOIx,VOIy,VOIz]=ind2sub(size(VOI),single(find(VOI)));
else
    voirect = getrect;
    recty = intersect(recty,fliplr(VOIdim(2) - round(voirect(1):(voirect(1)+voirect(3))))); % y is flipped in this view
    rectz = intersect(rectz,round(voirect(2):(voirect(2)+voirect(4))));

    [VOIx,VOIy,VOIz] = meshgrid(rectx, recty, rectz);
end
VOIvoxels = [VOIx(:) VOIy(:) VOIz(:)]; clear VOIx VOIy VOIz;

set(st.vols{1}.ax{3}.ax,axes_properties,restore_axes_values);
delete(texthandle);

cursorpos = spm_orthviews('Pos');
spm_orthviews('rmblobs',1)
if ndefaults.voi.reposition
    blobcenter = double(nut_coordtfm(mean(VOIvoxels),translation_tfm));
    spm_orthviews('Reposition',blobcenter);
end
spm_orthviews('addcolouredblobs',1,VOIvoxels',ones(size(VOIvoxels,1),1),translation_tfm,[.3 1 .2]);
st.vols{1}.blobs{1}.max=4;  % force desired color scale
st.vols{1}.blobs{1}.min=0;
spm_orthviews('redraw');

%%% translate VOIvoxels so that origin is at (0,0,0)mm
VOIvoxels = nut_coordtfm(VOIvoxels,translation_tfm);

VOIvoxels = nut_mri2meg(VOIvoxels);  % save in MEG head coords
% nuts.VOIvoxels = VOIvoxels;

if asksave
    [filename, pathname] = uiputfile('*.mat', 'Save VOI...(OPTIONAL)','voxels.mat');  %optional as long as saving VOIvoxels to nuts
    if ~(isequal(filename,0) | isequal(pathname,0))
        save(fullfile(pathname, filename),'VOIvoxels');
    end
end

% remove any lead field since we've obviously made it invalid by loading a new VOI
if(isfield(nuts,'Lp'))
    nuts = rmfield(nuts,'Lp');
end

%nut_enabler;
