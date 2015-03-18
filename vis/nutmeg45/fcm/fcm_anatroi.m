function [roi,beam] = fcm_anatroi(label,refimg,reftxt)
% Usages:
%  fcm_anatroi -list roidef
%     lists all ROI labels defined in ROI definition ROIDEF
%
%  roi = fcm_anatroi(label,roidef)
%     finds MNI voxel coordinates with selected LABEL in ROI definition ROIDEF
%
% roidef        The template to use:
%               'aal' (default) Automatic Anatomical Labelling
%               'mni' Montreal Neurological Institute
%               'hmat' Human Motor Area Template
%               or path to custom template image which must be in Analyze
%               format.

global nuts 

if nargin<2 || isempty(refimg)
    refimg='aal';
end

if strcmpi(label,'gui')
    if isempty(findobj('tag','nutmegfig'))
        if isempty(nuts)
            error('You must load session file first.')
        end
        nutmeg(nuts);
    end
    label = fcm_roiselect_gui;
    if isempty(label), roi=[]; return, end
    refimg = 'aal'; % other templates currently not supported for GUI selection.
end

if ~iscell(label), label={label}; end

switch lower(refimg)
case 'aal'
    refimg=which('ROI_MNI_V4.nii');
    if isempty(refimg), error('AAL toolbox must be in matlab path.'), end
    [pa,fi,ext]=fileparts(refimg);
    reftxt=fullfile(pa,[fi '.txt']);
    [l1,l2,ln]=textread(reftxt,'%s %s %d');
    
    if strcmp('-list',label{1})
        roi = [l1 l2];
        return
    end
    
    V = spm_vol(refimg);
    Y = spm_read_vols(V);
    voxelsize=diag(V.mat);
    voxelsize=abs(voxelsize(1:3));
    
    [vect,vb] = nut_vol2vector(Y);
    refvox = nut_coordtfm(vb,V.mat);
    
    f=[];
    for k=1:length(label)
        v = (strcmpi(label{k},l2) | strcmpi(label{k},l1));
        if ~any(v), error('label %s not found',label{k}), end
        label{k} = l2{v};  % make sure it's the long label
        f = [f; find(vect == ln(v))];
    end
    f=unique(f);
    
case 'mni'
    load('ihb_DataBase.cdb','-mat'); % load MNI labels
    
    if strcmp('-list',label{1})
        roi = MNIdb.cNames{5};
        return
    end
    
    [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
    refvox = [meshx(:) meshy(:) meshz(:)];
    voxelsize = [MNIdb.voxX MNIdb.voxY MNIdb.voxZ];
    
    f=[];
    for k=1:length(label)
        v = find(strcmpi(label{k},MNIdb.cNames{5}));
        if isempty(v), error('label %s not found',label{k}), end    
        f = [f; find(MNIdb.data(:,5) == v)];
    end
    f=unique(f);
            
case 'hmat'
    if isempty(which('BrikLoad.m')), error('AFNI Matlab Library is needed in matlab path to read BRIK file.'), end
    refimg=which('HMAT+tlrc.BRIK');
    if isempty(refimg), error('HMAT template must be in matlab path.'), end
    
    list={'M1 R';'M1 L';'S1 R';'S1 L';'SMA R';'SMA L';'pre-SMA R';'pre-SMA L'; ...
          'PMd R';'PMd L';'PMv R';'PMv L'};
    
    if strcmp('-list',label{1})
        roi = list;
        return
    end
    
    [~,Y,hdr] = BrikLoad(refimg);
    [vect,vb] = nut_vol2vector(Y);
    voxelsize = hdr.DELTA;
    mat = [-voxelsize(1) 0 0 -hdr.ORIGIN(1);0 -voxelsize(2) 0 -hdr.ORIGIN(2);0 0 voxelsize(3) hdr.ORIGIN(3);0 0 0 1];
    refvox = nut_coordtfm(vb,mat)-1;
    
    f=[];
    for k=1:length(label)
        v = find(strcmpi(label{k},list));
        if isempty(v), error('label %s not found',label{k}), end
        f = [f; find(vect == v)];
    end
    f=unique(f);

otherwise
    V = spm_vol(refimg);
    Y = spm_read_vols(V);
    
    [vect,vb] = nut_vol2vector(Y);
    refvox = nut_coordtfm(vb,V.mat);
    voxelsize=diag(V.mat);
    voxelsize=abs(voxelsize(1:3));    

    if ischar(label{1})
        [l1,l2,ln]=textread(reftxt,'%s %s %d');
        if strcmp('-list',label{1})
            roi = [l1 l2];
            return
        end
        
        f=[];
        for k=1:length(label)
            v = find(strcmpi(label{k},l2) | strcmpi(label{k},l1));
            if ~any(v), error('label not found'), end        
            f = [f; find(vect == ln(v))];
        end
        f=unique(f);
    else
        if strcmp('-list',label{1})
            roi = unique(vect);
            return
        end
        f=[];
        for k=1:length(label)
            f = [f; find(vect == label{k})];
        end
        f=unique(f);
    end
    
end

roi.MNIvoxels = refvox(f,:);
roi.label = label;

if isfield(nuts,'coreg') && isfield(nuts.coreg,'norm_mripath')
    fcm_roiidx(roi);
end

if nargout>1
    beam=struct('s',{{ones(size(roi.MNIvoxels,1),1)}},'voxels',roi.MNIvoxels,'voxelsize',voxelsize, ...
        'timepts',1,'timewindow',[0 2],'bands',[0 1],'srate',1);
    beam.coreg=struct('mripath',which('avg152T1.nii'),'orientation',1,'meg2mri_tfm',eye(4));
    
    save s_beam_roitemplate beam
    tv s_beam_roitemplate
end


