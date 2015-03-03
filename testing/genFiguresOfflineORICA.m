%% visualize ORICA results
EEG_orica = pop_loadset('20150115_Experiment_rmBadCH_AMICA.set');
load 20150203_result_rmBadCH_AMICA.mat;

% manually defined AMICA model regions
model_time =  [300, 700, 1100, 1500].*EEG_orica.srate;
time = cell2mat({results.time});
idx_old = 1;
for it = 1:length(model_time)
    idx = find(time>model_time(it),1);
    if ~isempty(idx) && idx ~= idx_old
        EEG_orica.icasphere = results(it).icasphere;
        EEG_orica.icaweights = results(it).icaweights;
        EEG_orica.icawinv = inv(EEG_orica.icaweights*EEG_orica.icasphere);
        
%         figure, pop_topoplot(EEG_orica,0,[1:30],sprintf('rmBadCH ORICA %d sec',time(idx)));
%         eval(sprintf('export_fig 20150203_rmBadCH_ORICA_comp_sess%d_1 -png -transparent',it));
%         close
%         
%         figure, pop_topoplot(EEG_orica,0,[31:58],sprintf('rmBadCH ORICA %d sec',time(idx)));
%         eval(sprintf('export_fig 20150203_rmBadCH_ORICA_comp_sess%d_2 -png -transparent',it));
%         close
        
        idx_old = idx;
    end
end


%% IC x Time x Correlation Plot - compare to Infomax
EEG_runica = pop_loadset('20150115_Experiment_raw_icainfo.set');
EEG_orica = pop_loadset('20150115_Experiment_raw_icainfo.set');
load 20150202_result_decayW8B8.mat;

nChs = length(results(1).icaweights);
nPts = length(results);
winv_runica = EEG_runica.icawinv_true;

% initialization
downFactor = 1;
sample = [1:downFactor:(nPts-1), nPts];
numCalc = length(sample);
corrComp = zeros(nChs,numCalc);
numPerm = zeros(1,numCalc);
time = zeros(1,numCalc);

hbar = waitbarTime(0,'Good results need patience...');
for it = 1:numCalc
    weights = results(sample(it)).icaweights * results(sample(it)).icasphere;
    winv = inv(weights);

    % compute correlation
    [correlation, ind1, ind2] = matcorr(winv_runica', winv', 0, 0); % ind1 and ind2 are ranked according to correlation
    [corrPerm,idx_corr] = sortrows([ind1, correlation],1);
    corrComp(:,it) = abs(corrPerm(:,2));

    % compute number of IC permutation compared to previous best match
    idx = sortrows([ind1, ind2],1);   % sort ind_orica based on ind_truth
    if it == 1
        numPerm(it) = sum(ind1 ~= ind2);
    else
        numPerm(it) = sum(bestmatch ~= idx(:,2));
    end
    bestmatch = idx(:,2);
    time(it) = results(sample(it)).time;
    waitbarTime(it/numCalc,hbar);   
end
close(hbar);

figure();
surf(time/EEG_orica.srate/60, 1:nChs, corrComp, 'EdgeColor','none');
xlabel('Time (min)'); ylabel('Component ID'); zlabel('corr coef'); %zlim([0 1]); 
ylim([1 nChs]); xlim([0 25.5]); view(2); 
hc = colorbar; title(hc,'Correlation','FontSize', 14, 'LineWidth',2,'Position', [0, -0.12]);

% plot correlation curve of selected component
selectComp = [1, 7, 12];
downfactor = 1; downIdx = 1:downfactor:length(time);
figure(); hold on
plot(time(downIdx)/EEG_runica.srate/60,corrComp(selectComp(1),downIdx),'g','LineWidth',2);
plot(time(downIdx)/EEG_runica.srate/60,corrComp(selectComp(3),downIdx),'b','LineWidth',2);
plot(time(downIdx)/EEG_runica.srate/60,corrComp(selectComp(2),downIdx),'r','LineWidth',2);
xlabel('Time (min)'); ylabel('Correlation coefficients'); ylim([0 1]);
legend('prefrontal','occipital','fronto-central');

% plot ground truth
selectComp = 1:30;
idx = sortrows([ind1, ind2],1);   % sort ind_orica based on ind_truth
winvPerm = winv(:,idx(:,2)) * diag(sign(correlation(idx_corr)));
EEG_orica.icawinv = winvPerm;
EEG_orica.icaweights = inv(winvPerm);
EEG_orica.icasphere = eye(nChs);
pop_topoplot(EEG_orica,0,selectComp);

% plot number of order-swtiching
figure;
plot(time/EEG_runica.srate/60, numPerm);
xlabel('Time (min)'); ylabel('Number of order switching');

% prettify_plot;
% export_fig 2015020_decayW8B8_ORICA-Infomax_corrmap -png -transparent;
% export_fig CorrMap_runica -png -transparent;
% export_fig CorrComp_smooth_runica -png -transparent;
% export_fig NumPerm_runica -png -transparent;


%% IC x Time x Correlation Plot - compare to AMICA
% specifically for 20150115_Experiment dataset with the AMICA model
EEG_AMICA = pop_loadset('filename','20150115_Experiment_AMICA.set','filepath','D:\\Matlab Coding\\VisEEG\\data\\');
load 20150203_result_rmBadCH_AMICA.mat;

% visualize AMIRCA model probability
models2plot = [1:5]; smoothed = 1; smoothlength = 2; store = 0;
EEG_AMICA = pop_modprobplot(EEG_AMICA,models2plot,smoothed,smoothlength,store);
% visualize AMIRCA model components
model = 1; components = [1:32];
pop_topohistplot(EEG_AMICA,0,components,'Luca_exp1, Model 1',0,0,'electrodes','on','showhist',0,'use_block',1,'model',model);

% store AMICA results
nModels = EEG_AMICA.etc.amica.num_models;
nChs = length(EEG_AMICA.icaweights);

% manually defined AMICA model regions
model_order = [2,   5,      3,      5,      4,      5,      1];
model_time =  [0,   330,    357,    717,    747,    1052,   1078].*EEG_AMICA.srate;

% initialization
nPts = length(results);
downFactor = 1; % downsampling
sample = [1:downFactor:(nPts-1), nPts]; 
numCalc = length(sample);
corrComp = zeros(nChs,numCalc);
numPerm = zeros(1,numCalc);
time = zeros(1,numCalc);

hbar = waitbarTime(0,'Good results need patience...');
for it = 1:numCalc
    time(it) = results(sample(it)).time;
    
    % assign ORICA model
    weights = results(sample(it)).icaweights * results(sample(it)).icasphere;
    winv = inv(weights);
    winv = winv(1:nChs,1:nChs); % match number of sources as AMICA model / need to fix this step

    % assign AMICA model
    model_idx = find(model_time<time(it));
    winv_amica = EEG_AMICA.etc.amica.A(:,:,model_order(model_idx(end)));
    
    % compute correlation
    [correlation, ind1, ind2] = matcorr(winv_amica', winv', 0, 0); % ind1 and ind2 are ranked according to correlation
    [corrPerm,idx_corr] = sortrows([ind1, correlation],1);
    corrComp(:,it) = abs(corrPerm(:,2));

    % compute number of IC permutation compared to previous best match
    idx = sortrows([ind1, ind2],1);   % sort ind_orica based on ind_truth
    if it == 1
        numPerm(it) = sum(ind1 ~= ind2);
    else
        numPerm(it) = sum(bestmatch ~= idx(:,2));
    end
    bestmatch = idx(:,2);
        
    waitbarTime(it/numCalc,hbar);   
end
close(hbar);

figure();
surf(time/EEG_AMICA.srate/60, 1:nChs, corrComp, 'EdgeColor','none');
xlabel('Time (min)'); ylabel('Component ID'); zlabel('corr coef'); %zlim([0 1]); 
ylim([1 nChs]); xlim([0 25.5]); view(2); 
hc = colorbar; title(hc,'Correlation','FontSize', 14, 'LineWidth',2,'Position', [0, -0.12]);

% export_fig 20150203_rmBadCH_ORICA-AMICA_corrmap -png -transparent;

%% Model fitness - PI x time for ORICA compared to each AMICA model
EEG_AMICA = pop_loadset('filename','20150115_Experiment_AMICA.set','filepath','D:\\Matlab Coding\\VisEEG\\data\\');
load 20150203_result_rmBadCH_AMICA.mat;

% store AMICA results
nModels = EEG_AMICA.etc.amica.num_models;
nChs = length(EEG_AMICA.icaweights);

% initialization
nPts = length(results);
downFactor = 1; % downsampling
sample = [1:downFactor:(nPts-1), nPts]; 
numCalc = length(sample);
modelPI = zeros(nModels,numCalc);
time = zeros(1,numCalc);

hbar = waitbarTime(0,'Good results need patience...');
for it = 1:numCalc
    time(it) = results(sample(it)).time;
    
    % assign ORICA model
    weights = results(sample(it)).icaweights * results(sample(it)).icasphere;
%     weights = weights(1:nChs,1:nChs); % wrong way to select components

    % compute model fitness
    for model_idx = 1:nModels
        H = weights * EEG_AMICA.etc.amica.A(:,:,model_idx);
        C = H.^2;
        modelPI(model_idx,it) = (nChs-sum(max(C,[],1)./sum(C,1))/2-sum(max(C,[],2)./sum(C,2))/2)/(nChs-1);
    end
    
    waitbarTime(it/numCalc,hbar);   
end
close(hbar);

figure();
plot(time/EEG_AMICA.srate/60,modelPI');
xlabel('Time (min)'); ylabel('Performance Index');
legend('Model 1','Model 2','Model 3','Model 4','Model 5');

% export_fig 20150203_rmBadCH_ORICA-AMICA_PI -png -transparent;

modelProb = zeros(nModels,numCalc);
for it = 1:numCalc
    for model_idx = 1:nModels
        modelProb(model_idx,it) = modelPI(model_idx,it) / sum(modelPI(:,it));
        modelProb(model_idx,it) = 1 - modelProb(model_idx,it);
    end
end

figure();
plot(time/EEG_AMICA.srate/60,modelProb');
xlabel('Time (min)'); ylabel('Model Probability');
legend('Model 1','Model 2','Model 3','Model 4','Model 5');

% export_fig 20150203_rmBadCH_ORICA-AMICA_PIProb -png -transparent;
