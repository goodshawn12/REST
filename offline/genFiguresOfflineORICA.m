%% visualize ORICA results
EEG_orica = pop_loadset('SIM_NSTAT_3sess_16ch_3min.set'); % 20150115_Experiment_rmBadCH_AMICA
% load sim_3sess_16ch_3min.mat; % 20150203_result_rmBadCH_AMICA.mat;

time = cell2mat({results.time});
lambdaMSP = cell2mat({results.lambda});
mir = cell2mat({results.mir});
figure,plot(time/EEG_orica.srate/60,10*log10(lambdaMSP)); ylabel('lambda');
figure,plot(time/EEG_orica.srate/60,mir); ylabel('mir');

if isfield(results,'normRn'),   normRnSMP = cell2mat({results.normRn}); figure,plot(time/EEG_orica.srate/60,normRnSMP); ylabel('normRn'); end
if isfield(results,'normRnW'), normRnWSMP = cell2mat({results.normRnW}); figure,plot(time/EEG_orica.srate/60,normRnWSMP); ylabel('normRnW'); end
if isfield(results,'lambdaGrad'), lambdaGradSMP = cell2mat({results.lambdaGrad}); figure,plot(time/EEG_orica.srate/60,lambdaGradSMP); ylabel('lambdaGrad'); end
if isfield(results,'normCovyy'), normCovyySMP = cell2mat({results.normCovyy}); figure,plot(time/EEG_orica.srate/60,normCovyySMP); ylabel('normCovyy'); end
if isfield(results,'normCovfy'), normCovfySMP = cell2mat({results.normCovfy}); figure,plot(time/EEG_orica.srate/60,normCovfySMP); ylabel('normCovfy'); end


%% Evaluate ORICA performance on simulated data
printfile = 0;
plotIC = 0;
% manually defined session boundary
model_time =  [180, 360, 540].*EEG_orica.srate; % model_time =  [300, 700, 1100, 1500].*EEG_orica.srate;

nsrc = EEG_orica.nbchan;
model = 1; 
LFM = EEG_orica.etc.LFM{model};
EEG_orica.icasphere = eye(nsrc);
EEG_truth = EEG_orica;

% Fig 1. Component maps of ORICA and ground truth for each session
convergence = zeros(1,length(results));
for it = 1:length(results)
    if time(it)>model_time(model)
        % permute icaweights based on corr, correlation coefficients.
        orica_weights = results(it-1).icaweights*results(it-1).icasphere;
        orica_icawinv = pinv(orica_weights);
        
        [correlation, ind_truth, ind_orica] = matcorr(LFM', orica_icawinv', 0, 0);
        idx_sort = sortrows([ind_truth, ind_orica],1);   % sort ind_orica based on ind_truth
        [~,idx_corr] = sortrows([ind_truth, correlation],1);
        EEG_orica.icawinv = orica_icawinv(:,idx_sort(:,2)) * diag(sign(correlation(idx_corr)));
        EEG_orica.icaweights = pinv(EEG_orica.icawinv);

        if plotIC, pop_topoplot(EEG_orica,0,1:nsrc,sprintf('ORICA %d min',time(it-1)/EEG_orica.srate/60),[2,nsrc/2]); end
        if printfile, eval(sprintf('export_fig sim_16ch_3min_orica_sess%d_turbo10 -png -transparent',model)); close; end
        
        EEG_truth.icawinv = LFM;
        EEG_truth.icaweights = pinv(LFM);
        if plotIC, pop_topoplot(EEG_truth,0,1:nsrc,sprintf('Ground Truth Sess %d',model),[2,nsrc/2]); end
        if printfile, eval(sprintf('export_fig sim_16ch_3min_truth_sess%d_turbo10 -png -transparent',model)); close; end

        if model < 3
            model = model+1;
            LFM = EEG_orica.etc.LFM{model};
        end
    end
    H = results(it).icaweights * results(it).icasphere * LFM;
    C = H.^2;
    convergence(it) = (nsrc-sum(max(C,[],1)./sum(C,1))/2-sum(max(C,[],2)./sum(C,2))/2)/(nsrc-1);
end

% Fig 2. Performance Index over time
figure, plot(time/EEG_orica.srate/60,10*log10(convergence));
ylabel('Performance Index (dB)'); xlabel('Time (min)'); 
if printfile, export_fig sim_16ch_3min_PI_turbo -png -transparent; end

t = 1:length(normRn);
figure,
subplot(4,2,1); plot(t/EEG_orica.srate/60,10*log10(mean(lambda,1))); ylabel('lambda');
subplot(4,2,2); plot(time/EEG_orica.srate/60,10*log10(convergence)); ylabel('PI');
subplot(4,2,3); plot(time/EEG_orica.srate/60,mir); ylabel('mir');
subplot(4,2,4);plot(t/EEG_orica.srate/60,lambdaGrad); ylabel('lambdaGrad'); ylim([min(lambdaGrad(1*60*128:end)), max(lambdaGrad(1*60*128:end))]);
subplot(4,2,5);plot(t/EEG_orica.srate/60,normRn); ylabel('normRn'); ylim([min(normRn(1*60*128:end)), max(normRn(1*60*128:end))]);
subplot(4,2,6);plot(t/EEG_orica.srate/60,normRnW); ylabel('normRnW'); ylim([min(normRnW(1*60*128:end)), max(normRnW(1*60*128:end))]);
subplot(4,2,7);plot(t/EEG_orica.srate/60,normCovyy); ylabel('normCovyy'); ylim([min(normCovyy(1*60*128:end)), max(normCovyy(1*60*128:end))]);
subplot(4,2,8);plot(t/EEG_orica.srate/60,normCovfy); ylabel('normCovfy'); ylim([min(normCovfy(1*60*128:end)), max(normCovfy(1*60*128:end))]); xlabel('time(min)');
% export_fig sim_adapt_Murata_success -png -transparent;



%% Fig 3. Correlation between ORICA-decomposed ICs and the ground truth over time
nPts = length(results);
downFactor = 1; % downsampling
sample = [1:downFactor:(nPts-1), nPts]; 
numCalc = length(sample);
corrComp = zeros(nsrc,numCalc);
numPerm = zeros(1,numCalc);
timeSample = zeros(1,numCalc);

model_Sess = [0,23040,46080]; model = 0;
hbar = waitbarTime(0,'Good results need patience...');
for it = 1:numCalc
    timeSample(it) = results(sample(it)).time;
    
    % assign ORICA model
    weights = results(sample(it)).icaweights * results(sample(it)).icasphere;
    winv = inv(weights);

    % assign ground truth
    model_new = find(timeSample(it)>model_Sess,1,'last');
    if model ~= model_new
        winv_truth = EEG_orica.etc.LFM{model_new};
        model = model_new;
    end
    
    % compute correlation
    [correlation, ind1, ind2] = matcorr(winv_truth', winv', 0, 0); % ind1 and ind2 are ranked according to correlation
    [corrPerm,idx_corr] = sortrows([ind1, correlation],1);
    corrComp(:,it) = abs(corrPerm(:,2));

    % compute number of IC permutation compared to previous best match
    idx = sortrows([ind1, ind2],1);
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
surf(timeSample/EEG_orica.srate/60, 1:nsrc, corrComp, 'EdgeColor','none');
xlabel('Time (min)'); ylabel('Component ID'); zlabel('corr coef'); %zlim([0 1]); 
ylim([1 nsrc]); xlim([0 9]); view(2); 
hc = colorbar; title(hc,'Correlation','FontSize', 14, 'LineWidth',2,'Position', [0, -0.12]);
if printfile, export_fig sim_16ch_3min_ORICA_corrmap_len250 -png -transparent; end;



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
