function [nuts, Lp, voxels]= nut_obtain_lead_field(nuts,voxelsize,lfcomp,forwfile,voxelfile);
% function [nuts, Lp, voxels]= nut_obtain_lead_field(nuts,voxelsize,lfcomp,forwfile,voxelfile);


if isfield(nuts,'meg')
    if ~isfield(nuts.meg,'system')
        % legacy support for old session.mat files with data loaded before
        % this was added, so we assume that the traditional
        % nut_compute_lead_field is still the one to call, without figuring
        % out which manufacturer it was (as CTF, BTi, and KIT all were okay
        % with it)
        nuts.meg.system='unknown';
    end
else
    error('you need to load some data first');
end

% NUTEEG mod: compute lead potentials
if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
    [Lp,voxels]=nut_compute_lead_potential(voxelsize);
    nuts.voxelsize=[voxelsize voxelsize voxelsize];
    return;
end


switch nuts.meg.system
    case {'CTF','CTF_orig','BTi','KIT','4D_BTi','yokogawa160','unknown'}
        % FIXME make above labels match Fieldtrip labels ft_senstype(data)
        if ~exist('lfcomp','var')
            lfcomp = 3;
        end
        if nuts.meg.grad_order & nuts.meg.keepref_ts
            chanmixMtx=nuts.meg.chanmixMtx{2}*nuts.meg.chanmixMtx{1};
        else
            chanmixMtx=nuts.meg.chanmixMtx{1};
        end
        [Lp,voxels]=nut_compute_lead_field(voxelsize,lfcomp,chanmixMtx);
        nuts.Lp=Lp;
        nuts.voxelsize=[voxelsize voxelsize voxelsize];
    case 'Neuromag'
        disp('this option allows to load leadfield from elsewhere')
        if ~exist('forwfile','var')
            [forw_filename, forw_path]=uigetfile('*.fif;*.bin','Select *.fif (MNE) or *.bin (Brainstorm) file with lead field included');
            if isequal(forw_filename,0)|isequal(forw_path,0)
                return;
            end
            forwfile=[forw_path forw_filename];
        end
        [crap,morecrap,ext] = fileparts(forwfile);
        switch(ext)
            case '.fif'  % call MNE code
                nuts=nut_mne2nuts(nuts,forwfile);
            case '.bin'  % call Brainstorm code
                if ~exist('voxelfile','var')
                    [vox_filename, vox_path]=uigetfile('*.mat*','Select file with voxel locations');
                    if isequal(vox_filename,0)|isequal(vox_path,0)
                        return;
                    end
                end
                voxelfile=[vox_path vox_filename];
                nuts=nut_bst2nuts(nuts,forwfile,voxelfile);
            otherwise
                error('Are you loading a lead field created outside of MNE or Brainstorm?  Please contact developers');
        end
        Lp=nuts.Lp;
        voxels=nuts.voxels;
        nuts.voxelsize=[voxelsize voxelsize voxelsize];
    case 'BrainProducts'
        sens.type='eeg';
        nuts.voxelsize=[voxelsize voxelsize voxelsize];
        [nuts,Lp]=nut_fieldtrip_compute_leadfield(nuts,sens);
        nuts.Lp=Lp;
        voxels=nuts.voxels;
    case 'Neuroscan'
        [Lp,voxels]=nut_compute_lead_potential(voxelsize);
        nuts.Lp=Lp;
        nuts.voxelsize=[voxelsize voxelsize voxelsize];
    case 'matfile'
        error('change nuts.meg.system to the manufacturer type, or load in your lead field explicitly');
    otherwise
        error('datatypes here should match ones in nut_importmeg...(ask JZ?)')
end

