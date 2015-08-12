%{
make_headModel    Wrapper function to generate MoBILAB headModel object to 
                  use with REST.

*Requires openMEEG*
    
HMOBJ = makeHeadModel(NAME,CHANLOCS,OUTPUT_DIR,SURFACES_FILE,CONDUCTIVITY)

Inputs:
    NAME:           Name to use for output files (defualt: none)
    CHANLOCS:       EEGLAB format channel locations structure (defualt: none)
    NORMAL2SURFACE: True if the model constrains sources in the cortex to
                    be normal to the cortical surface
                    False if sources are free to point in any direction
                    (default: false)
    OUTPUT_DIR:     Directory to which output files will be saved 
                    (default: current working directory)
    SURFACES_FILE:  Matfile containing information about the head model to be made
                     - BEM model in variable "surfData" as 1x3 structure with
                       fields: vertices: vertices of the mesh (nPoints x 3)
                               faces: triangular faces of the mesh (nFaces x 3)
                       Elements shoud go from out to in
                       (e.g., skin - skull - cortex)
                     - Fiducials in variable "fiducials" as a structure
                       with fields: nason, lpa, rpa, vertex, inion. All are
                       1x3 numerical arrays giving the position in space of
                       that fiducial
                     - Atlas (optional) in variable "atlas" as a structure
                       with fields: colorTable: integer index of regions
                                    label: cellstr array of region names
    CONDUCTIVITY:   Conductivity (S/m) values to use with the BEM mesh.
                    Values match the order of surfData in the surfaces_file
                                    
Outputs:
    HMOBJ:          headModel object generated. Also stored in a file.

*Requires openMEEG*
%}
function hmObj = make_headModel(name,chanlocs,normal2surface,output_dir,surfaces_file,conductivity)
%% check inputs
if ~ischar(name)
    error('make_headModel: name must be a string'); end
if ~exist('chanlocs','var') || isempty(chanlocs)
    error('make_headModel: no channel locations provided'); end
if ~exist('output_dir','var') || isempty(output_dir)
    disp('make_headModel: no output directory provided, using current directory')
    output_dir = './';
end
if ~exist('normal2surface','var') || isempty(normal2surface)
    normal2surface = false; end
if ~exist('conductivity','var') || isempty(conductivity)
    conductivity = [0.33 0.022 0.33]; % brain and scalp = 0.33 S/m, skull = 0.022 S/m; these conductivies were taken from
                                  % Valdes-Hernandez et al., 2009, Oostendrop TF, 2000; Wendel and Malmivuo, 2006
end


%% load meshes and other data
if ~exist('surfaces_file','var') || isempty(surfaces_file)
    % use Colin27 head
    disp('Using Colin27 head mesh')
    mobilabPath = fileparts(which('headModel'));
    templateFile = [mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat'];
    template = load([mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat']);
    surfFile = [name '_Colin27_4825.mat'];
    surfData = template.surfData;
else
    % use provided head
    disp('Using provided head mesh')
    templateFile = [mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat'];
    template = load([mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat']);
    surfFile = [name '_Colin27_4825.mat'];
    surfData = template.surfData;
end


%% display template head model: Colin27 + aal atlas + emotiv montage
dipfitdefs;
vol.bnd.pnt = surfData(3).vertices;
vol.bnd.tri = surfData(3).faces;
index = num2str(randi(1e5));
save(['temp' index],'vol')
% rmpath(genpath([mobilabPath '/dependency/iso2mesh']))
electrodes = coregister(chanlocs,template_models(2).chanfile,'mesh',['temp' index]);
% addpath(genpath([mobilabPath '/dependency/iso2mesh']))
delete(['temp' index '.mat'])

elec = electrodes.pnt;
label = {chanlocs.labels};
if isfield(template,'atlas')
    hmObj = headModel('surfaces',surfFile,'fiducials',template.fiducials,'channelSpace',elec,'label',label);
else  
    hmObj = headModel('surfaces',surfFile,'atlas',template.atlas,'fiducials',template.fiducials,'channelSpace',elec,'label',label);
end


%% warping channel space to template
aff = hmObj.warpChannelSpace2Template(templateFile,surfFile,'affine');


%% solving the forward problem with OpenMEEG

dir_original = cd;
cd(output_dir);

hmObj.computeLeadFieldBEM(conductivity,normal2surface);

cd(dir_original);


%% solving the inverse problem with sLORETA
% K: lead field matrix
% L: Laplacian operator
% rmIndices: indices to be removed (the Thalamus)
% surfData(3).vertices: source space



%--
[U,S,V] = svd(K/L,'econ');
Ut = U';
s2 = diag(S).^2;
iLV = L\V;



eval([name '_HeadModel = hmObj;']);
save([output_dir filesep name '_HeadModel'],[name 'HeadModel'])
% save headModel
hmObj.saveToFile([fileparts(output_dir) filesep name])

% save LFM and Laplacian for cropped cortex
[~,K,L,rmIndices] = getSourceSpace4PEB(hmObj);
hmInd = setdiff(1:handles.nVertices,rmIndices);
save([output_dir filesep name '_SSPEB'],'K','L','hmInd')
