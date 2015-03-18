for kk=1:5
    switch('mne')
        case 'brainstorm'
           % gain = lead field = forward model = forward field = green's matrix = .....
           % Lxyz is 3-component lead field from BrainStorm
           % need to use brainstorm/headmodel/read_gain.m to read
           % 3-component ***non-matfile*** brainstorm lead field (*.bin)
           % resave as matlab matrix 'Lxyz.mat'
           nuts = nut_bst2nuts(['raw_hc080251_090402_run' num2str(kk) '_tsss.fif'],'Lxyz','tess_cortex_concat_15000V');
        case 'mne'
           % but brainstorm's lead field didn't work for us despite our sincere best efforts in one dataset, so use MNE
           nuts = nut_mne2nuts(['raw_hc080251_090402_run' num2str(kk) '_tsss.fif'],'GrandAverageH_C-7-fwd.fif');
    end
    nuts = nut_neuromagepoch(nuts,11,[-.25 .7]);  % 11 is event value; [-.25 .7] is window around marker to segment in seconds
    nuts.meg.data=nuts.meg.data2;
    nutrun{kk}=nuts;
    clear nuts
end

nuts = nutrun{1};
nuts.meg = rmfield(nuts.meg,'data2');
nuts.Lp = double(nuts.Lp);

% concatenate all runs into one nutmeg session file
nuts.meg.data = cat(3,[nutrun{:}.meg.data2]);

save nutsession-event11 nuts;

% nutrun{2}.meg.data2,nutrun{3}.meg.data2,nutrun{4}.meg.data2,nutrun{5}.meg.data2)