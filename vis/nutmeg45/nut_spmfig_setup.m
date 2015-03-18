function nut_spmfig_setup
% adds custom features to SPM MRI display window
% e.g., MEG coords, MNI coords, activation intensity

global st nuts beam rivets ndefaults

if(~isempty(beam))
    voxelsize = beam.voxelsize; % this determines scroll bar increments...
    coreg = beam.coreg;
else
    voxelsize = [5 5 5];
    coreg = nuts.coreg;
    rivets.sliderenable = 'on';
end

fg = spm_figure('GetWin','Graphics');
WS = spm('winscale');

% reposition MRI intensity box
inlabel = findobj('String','Intensity:','HorizontalAlignment','center');

% get rid of some stuff, sweep other stuff under the rug
delete(findobj('String','Crosshair Position'));
delete(findobj('ToolTipString','move crosshairs to origin'));
delete(inlabel);
set(st.in,'Visible','off');

beaminlabel = uicontrol(fg,'Style','Text','Position',[60 350 120 020].*WS,'String','Activation Intensity:');
st.beamin = uicontrol(fg,'Style','edit', 'Position',[175 350  85 020].*WS,'String','');

uicontrol(fg,'Style','Text', 'Position',[75 255 35 020].*WS,'String','MNI:');
st.mnip = uicontrol(fg,'Style','edit', 'Position',[110 255 135 020].*WS,'String','','Callback','nut_image(''setposmni'')','ToolTipString','move crosshairs to MNI mm coordinates');
st.mnilabel=uicontrol('Style','text','BackgroundColor',[1 1 1],'Units','normalized','Position',[.6 .5 .3 .15],'FontSize',12);

uicontrol(fg,'Style','Text', 'Position',[75 315 35 020].*WS,'String','MEG:');
st.megp = uicontrol(fg,'Style','edit', 'Position',[110 315 135 020].*WS,'String',sprintf('%.1f %.1f %.1f',nut_mri2meg(spm_orthviews('pos')')),'Callback','nut_image(''setposmeg'')','ToolTipString','move crosshairs to MEG mm coordinates');

set(st.mp,'Callback','spm_image(''setposmm''); nut_image(''shopos'');');
set(st.vp,'Callback','spm_image(''setposvx''); nut_image(''shopos'');');

%add sliders for scrolling through MRI
if ndefaults.sliders
    ax_pos = get(st.vols{1}.ax{1}.ax,'Position');
    cor_pos = get(st.vols{1}.ax{2}.ax,'Position');
    sag_pos = get(st.vols{1}.ax{3}.ax,'Position');
    slidercallback1 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([st.vols{1}.dim(1)*get(gcbo,''Value'') pos(2) pos(3)])); spm_image(''setposvx''); nut_image(''shopos'')';
    slidercallback2 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) st.vols{1}.dim(2)*get(gcbo,''Value'') pos(3)])); spm_image(''setposvx''); nut_image(''shopos'')';
    slidercallback3 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) pos(2) st.vols{1}.dim(3)*get(gcbo,''Value'')])); spm_image(''setposvx''); nut_image(''shopos'')';
    
    mat = st.vols{1}.premul*st.vols{1}.mat;
    R=spm_imatrix(st.vols{1}.mat);
    R = spm_matrix([0 0 0 R(4:6)]);
    R = R(1:3,1:3);
    dim = st.vols{1}.dim*inv(R);
    dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
    
    st.slider{1}=uicontrol(fg,'Style','Slider','Callback',slidercallback1,'SliderStep',[voxelsize(1)/abs(dim_mm(1)) .1],'Units','normalized','Position',[ax_pos(1) ax_pos(2)+ax_pos(4) ax_pos(3) .015]);
    st.slider{2}=uicontrol(fg,'Style','Slider','Callback',slidercallback2,'SliderStep',[voxelsize(2)/dim_mm(2) .1],'Units','normalized','Position',[ax_pos(1)+ax_pos(3) ax_pos(2) .02 ax_pos(4)]);
    st.slider{3}=uicontrol(fg,'Style','Slider','Callback',slidercallback3,'SliderStep',[voxelsize(3)/dim_mm(3) .1],'Units','normalized','Position',[cor_pos(1)+cor_pos(3) cor_pos(2) .02 cor_pos(4)]);
    set([st.slider{:}],'Visible',rivets.sliderenable);
end

%                 set(st.slider{i},'Value',max((posmrimm(i)-st.bb(1,i))/dim(i),0));

% tell time series to refresh after clicking around in SPM volume
% rivets.ts_refresh_handle = @plot_ts;  %% pass handle since technically it's a private function

for i=1:3  % step through three orthogonal views
    SPM_axes_obj(i) = st.vols{1}.ax{i}.ax;
end

if(strcmp(spm('ver'),'SPM2'))
    if(isempty(findstr(get(SPM_axes_obj(1),'ButtonDownFcn'),'nut_image')))
        SPM_axes_ButtonDownFcn = [get(SPM_axes_obj(1),'ButtonDownFcn') 'nut_image(''shopos'');'];
        set(SPM_axes_obj,'ButtonDownFcn',SPM_axes_ButtonDownFcn);
    end
elseif(strcmp(spm('ver'),'SPM8b') || strcmp(spm('ver'),'SPM8'))
    nut_image_handle = @nut_image;
    tmp = get(SPM_axes_obj(1),'ButtonDownFcn');
    if(isa(tmp,'function_handle')) % this happens only with raw SPM8 ButtonDownFcn
        rivets.spm_refresh_handle = tmp;
    end
    SPM_axes_ButtonDownFcn = {@nut_image,'shopos'};
    set(SPM_axes_obj,'ButtonDownFcn',SPM_axes_ButtonDownFcn);
end

[dum,mrifile,dum]=fileparts(coreg.mripath);
if( isfield(coreg,'norm_mripath') || strcmp(mrifile(1),'w') || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )
%     load('ihb_DataBase.cdb','-mat'); % load MNI labels
%     rivets.MNIdb = MNIdb;
%     [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
%     rivets.MNIdb.coords = [meshx(:) meshy(:) meshz(:)];
    
    rivets.TalDB = load('talairachDB');
end
