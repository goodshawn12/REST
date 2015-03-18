%%
% Computes the lead potential.
% @param desiredvoxelsize size of each voxel.
% @param VOIvoxels the volume of interest.
% @return Lp lead potentials in V. [chan x vec x vox]
% @return voxels voxel coordinates.
% @author Daniel D.E. Wong
%%
function [Lp,voxels]=nut_compute_lead_potential(voxelsize)

global nuts

hWb = waitbar(0,'Please wait... computing lead potential...');

if ~isfield(nuts,'voxels') || (nuts.voxelsize(1)~=voxelsize(1))
    voxels=nut_make_voxels(voxelsize);

    answer=menu('Do you want to strip voxels?','No','Grey Matter','Inner Volume','Cortical Surface');
    if answer>1
        nuts.voxels=voxels;
        switch answer
            case 2
                nut_stripgreyvoxels;
            case 3
                nut_stripVoxels;
            case 4
                nut_stripcsvoxels;
        end
        voxels=nuts.voxels;
    end
else
    disp('Using existing nuts.voxels and nuts.voxelsize...')
    voxels = nuts.voxels;
end

Lp = zeros(size(nuts.meg.sensorCoord,1),3,size(voxels,1));
hdm = questdlg('Select head model','Compute Lead Potential','Spherical','Helsinki BEM','OpenMEEG BEM','Spherical');
if strcmp(hdm,'Spherical')
    if ~isfield(nuts.meg,'lsc'); error('Sphere center not set.'); end;
    radmeth = questdlg('Select radii calculation method:','Multisphere Computation','Standard Radii','Radii From Meshes','Standard Radii');
    
    if strcmp(radmeth,'Radii From Meshes');
        file = 0;
        while isequal(file,0); [file,path] = uigetfile('standard_vol.mat','Open Mesh File'); end;
        load(strcat(path,file));
        % Convert to MEG coordinates
        for i = 1:length(vol.bnd)
            vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
        end
        disp('Sphere radii will be computed from meshes')
    end
    
    tic

    lsc = nuts.meg.lsc;

    if(size(lsc,1)==1)
        disp('Single sphere')
        
        for ch = 1:size(nuts.meg.sensorCoord,1)
            if exist('vol','var')
                [lx1,ly1,lz1] = nut_bergsphleadp(voxels,nuts.meg.sensorCoord(ch,:),lsc,vol);
            else
                [lx1,ly1,lz1] = nut_bergsphleadp(voxels,nuts.meg.sensorCoord(ch,:),lsc);
            end
            Lp(ch,1,:) = lx1; Lp(ch,2,:) = ly1; Lp(ch,3,:) = lz1;
            waitbar(ch/size(nuts.meg.sensorCoord,1),hWb);
        end

    else    % Multispheres - need to calculate each channel one by one
        disp('Multispheres')
        
        % map lsc entries to channels...
        for i=1:size(nuts.meg.sensorCoord,1)
            indexmap(i)=strmatch(nuts.meg.sensor_labels{i},nuts.meg.lsc_sensor_labels,'exact');
        end

        for ch = 1:size(nuts.meg.sensorCoord,1)
            if exist('vol','var')
                [lx1,ly1,lz1] = nut_bergsphleadp(voxels,nuts.meg.sensorCoord(ch,:),lsc(indexmap(ch),:),vol);
            else
                [lx1,ly1,lz1] = nut_bergsphleadp(voxels,nuts.meg.sensorCoord(ch,:),lsc(indexmap(ch),:));
            end
            Lp(ch,1,:) = lx1; Lp(ch,2,:) = ly1; Lp(ch,3,:) = lz1;
            waitbar(ch/size(nuts.meg.sensorCoord,1),hWb);
        end

    end
    toc

elseif strcmp(hdm,'Helsinki BEM')    % Compute 3 layer BEM - Based on code provided by Matti Stenroos, 2008
    if(~exist('SetBEMPaths','file'))
        errordlg('The Helsinki BEM Library is required in the MATLAB path. It can be downloaded from http://peili.hut.fi/BEM')
        return
    end
    
    [file,path] = uigetfile('*_vol.mat','Open Head Surface Mesh File (required)');
    if isequal(file,0) || isequal(path,0)
        return
    end    
    load(strcat(path,file));
    
    % Convert to MEG coordinates
    for i = 1:length(vol.bnd)
        vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
    end
    
    [file,path] = uigetfile('*_bem.mat','Open BEM Matrices if available (cancel otherwise)');
    if ~isequal(file,0) & ~isequal(path,0)
        load(strcat(path,file));
    end
    
    tic

    if ~exist('meshes','var')
        % Prepare meshes, in [m]
        for i = 1:length(vol.bnd)
            meshes{i}=PrepareTriangleMesh(vol.bnd(i).vertices/1000,vol.bnd(i).faces(:,[1 3 2]));    %scalp->brain
        end

        %Conductivity information from the original model file [S/m]
        ci=vol.cond;
        co=[0 ci(1:end-1)];

        omega=DLMatrices_Linear(meshes);

        %The BEM matrix needed in the ISA approach: H�m�l�inen-Sarvas, IEEE TBME 36:165--171, Feb 1989;
        L0=TransferMatrix_Phi_Linear(meshes(end),ci(end),0,1,1,omega(end,end));
        %The standard BEM matrix (conductivity terms of Eq. (16) not multiplied inside the matrix)
        [L,inds,b,c]=TransferMatrix_Phi_Linear(meshes,ci,co,1,0,omega);

        [file,path] = uiputfile('_bem.mat','Save BEM matrices');  % BEM points are saved in MRI coordinates
        if ~isequal(file,0) & ~isequal(path,0)
            save(strcat(path,file),'meshes','ci','co','L0','L','inds','b','c');
        end
    end
    
    % Project electrodes to scalp triangle with closest centroid
    scalpmesh=PrepareTriangleMesh(vol.bnd(1).vertices,vol.bnd(1).faces(:,[1 3 2]));    %scalp->brain
    sensorCoord = zeros(size(nuts.meg.sensorCoord));
    senstrimap = zeros(size(nuts.meg.sensorCoord,1),1);
    sensorweights = zeros(size(nuts.meg.sensorCoord));
    for i = 1:size(nuts.meg.sensorCoord,1)
        senstrimap(i) = dsearchn(scalpmesh.mp,nuts.meg.sensorCoord(i,:));
        sensorCoord(i,:) = nuts.meg.sensorCoord(i,:) - scalpmesh.n(senstrimap(i),:).*dot(scalpmesh.n(senstrimap(i),:),nuts.meg.sensorCoord(i,:)-scalpmesh.mp(senstrimap(i),:))/sum(scalpmesh.n(senstrimap(i),:).^2);
        
        edge12 = scalpmesh.p(scalpmesh.e(senstrimap(i),2),:) - scalpmesh.p(scalpmesh.e(senstrimap(i),1),:);
        edge23 = scalpmesh.p(scalpmesh.e(senstrimap(i),3),:) - scalpmesh.p(scalpmesh.e(senstrimap(i),2),:);
        edge13 = scalpmesh.p(scalpmesh.e(senstrimap(i),3),:) - scalpmesh.p(scalpmesh.e(senstrimap(i),1),:);
        edgeS1 = sensorCoord(i,:)-scalpmesh.p(scalpmesh.e(senstrimap(i),1),:);
        edgeS2 = sensorCoord(i,:)-scalpmesh.p(scalpmesh.e(senstrimap(i),2),:);
        sensorweights(i,1) = norm(cross(edge23,edgeS2));
        sensorweights(i,2) = norm(cross(edge13,edgeS1));
        sensorweights(i,3) = norm(cross(edge12,edgeS1));
        sensorweights(i,:) = sensorweights(i,:)/sum(sensorweights(i,:));
    end
    
    if length(vol.bnd) == 3; disp('3 layer model detected. Will use isolated skull approach.'); end;
    
    for i = 1:size(voxels,1) % dsearchn(voxels,[20 -40 50])
        dpos = voxels(i,:)/1000;    % Dipole location in [m]
        for j = 1:3
            dmom = zeros(1,3);
            dmom(j) = 1e-9; % 1 nAm source
            
            %Calculate infinite medium potentials & solution with ISA
            for I=1:length(vol.bnd)
                phiinf{I}=DipolePotential_Infinite(meshes{I}.p,dpos,dmom);
                g{I}=phiinf{I}*b(I);%Eq. (16)
                if length(vol.bnd)==3; h{I}=g{I}/ci(end); end;  % ISA
                if length(vol.bnd)==3 & I==3                    % ISA
                    W0_O=L0*phiinf{end};
                    h{I}=h{I}-2/(ci(end)+co(end))*W0_O;
                end
            end
            hvec = [];
            if length(vol.bnd) == 3                             % ISA
                for k = 1:length(h); hvec = [hvec; h{k}]; end;
                hvec=co(end)*hvec;
                W=L*hvec;
                V=W;
                V(inds{end})=W(inds{end})+W0_O;
            else
                V=L*g{1};     %;g{2};g{3}
            end
            
%             global Vbem
%             Vbem = V(inds{1});
%             return;
            
            for k = 1:size(nuts.meg.sensorCoord,1)
                Lp(k,j,i) = sum(V(meshes{1}.e(senstrimap(k),:)).*sensorweights(k,:)');
            end
        end
        waitbar(i/size(voxels,1),hWb);
    end
    toc
elseif strcmp(hdm,'OpenMEEG BEM')
    if exist('ft_leadfield_openmeeg','file')
        [Lp,voxels] = nut_ftopenmeeg(voxelsize);
    else
        [Lp,voxels] = nut_om_compute_lfp(voxelsize);    % Native OpenMEEG pipeline
    end
end

waitbar(1,hWb);
delete(hWb);


%%
% Computes N nearest points to XI from X.
% @param X list of points.
% @param XI reference point.
% @param N number of nearest points to return.
% @return k indexes of nearest points.
%%
function k = dsearch2(X,XI,N)
k = zeros(N,1);
infvec = Inf*ones(size(XI));
for n = 1:N
    k(n) = dsearchn(X,XI);
    X((k(n)),:) = infvec;
end