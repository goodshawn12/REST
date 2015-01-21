% To build the head model from the files you have
T = load('../cognionicsHeadModel.mat')
hmObj = T.cognionicsHeadModel
hmObj.surfaces = '/home/lpiontonachini/Desktop/eeglab/plugins/mobilab/data/head_modelColin27_4825.mat';%mobilab.preferences.eeg.headModel
hmObj.computeLeadFieldBEM();

% Save the head model as an object
hmObj.saveToFile('headmodel_Luca_64.mat')

% Load the head model as an object
hmObj = headModel.loadFromFile('headmodel_Luca_64.mat');
hmObj.leadFieldFile = [hmObj.leadFieldFile '.mat']


[sourceSpace,K,L,rmIndices] = getSourceSpace4PEB(hmObj);

[U,S,V] = svd(K/L,'econ');
Ut = U';
s2 = diag(S).^2;
iLV = L\V;

load(hmObj.surfaces)
t = (0:512-1)/512;
x = cos(2*pi*10*t);
n = size(surfData(3).vertices,1);
ind = setdiff(1:n,rmIndices);


%%   This is your code
% simulating a Gaussian sources
I = strfind(hmObj.atlas.label,'Precuneus');
I = ~cellfun(@isempty,I);
roi = find(I);
hmObj.atlas.label(roi)
roi = hmObj.atlas.colorTable == roi(end);
roi = surfData(3).vertices(roi,:);
roi = mean(roi);
d = sqrt(sum((surfData(3).vertices(ind,:) - ones(length(ind),1)*roi).^2,2));
gSource = normpdf(d,0,10);
gSource = 0.001*gSource/max(gSource);
gSource = bsxfun(@times,gSource,x);
J = zeros(n,length(t));
Jtrue = J;
Jest = J;
Jtrue(ind,:) = gSource;

%-- Create a Jtrue as Nsources X 3
K_full = load(hmObj.leadFieldFile,'K');
K_full = K_full.K;
Jtrue = [Jtrue;Jtrue;Jtrue];
Jest  = 0*Jtrue;
Vtrue = K_full*Jtrue;

% repeat ind so it can index the other components  
I = false(n,3);
I(ind,:) = true;
ind = find(I(:));

hmObj.plotOnModel(Jtrue,Vtrue,'True source');

nlambdas = 100;
plotCSD = false;
Jest(ind,:) = inverseSolutionLoreta(Vtrue,K,L,nlambdas,plotCSD);
Vest = K*Jest(ind,:);
hmObj.plotOnModel(Jest,Vest,'Estimated source (Loreta)');

[Jest(ind,:),alpha,beta] = variationalDynLoreta(Ut,Vtrue,s2,iLV,L);
Vest = K*Jest(ind,:);
hmObj.plotOnModel(Jest,Vest,'Estimated source (varDynLoreta)');

[Jest(ind,:),alpha,beta] = dynamicLoreta(Ut,Vtrue,s2,iLV,L);
Vest = K*Jest(ind,:);
hmObj.plotOnModel(Jest,Vest,'Estimated source (varDynLoreta)');
%--

