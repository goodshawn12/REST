function labels = run_ICL(images, psds)

% get path information
pluginpath = fileparts(which('pop_iclabel'));

% activate matconvnet
run(fullfile(pluginpath, 'matconvnet', 'matlab', 'vl_setupnn'))

% load network
netStruct = load(fullfile(pluginpath, 'netICL.mat'));
net = dagnn.DagNN.loadobj(netStruct);
clear netStruct;

% inference
images = cat(4, images, -images, images(:, end:-1:1, :, :), -images(:, end:-1:1, :, :));
psds = repmat(psds, [1 1 1 4]);
try
    net.eval({
        'in_image', single(images), ...
        'in_psdmed', single(psds)
    });
catch e_run
    disp 'Failed to run ICLabel. Trying to compile MEX-files.'
    curr_path = pwd;
    cd(fileparts(which('vl_compilenn')));
    try
        vl_compile_nn
        cd(curr_path)
        disp(['MEX-files successfully compiled. Attempting to run ICLabel again. ' ...
            'Please consider emailing Luca Pion-Tonachini at lpionton@ucsd.edu to ' ...
            'share the compiled MEX-files. They will likely help other EEGLAB users ' ...
            'with similar computers as yourself.'])
        net.eval({
            'in_image', single(images), ...
            'in_psdmed', single(psds)
        });
    catch e_
        cd(curr_path)
        error(['MEX-file compilation failed. Further instructions on compiling ' ...
               'the MEX-files can be found at http://www.vlfeat.org/matconvnet/install/. ' ...
               'Further, you may contact Luca Pion-Tonachini at lpionton@ucsd.edu for help. ' ...
               'If you solve this issue without help, please consider emailing Luca as the ' ...
               'compiled files will likely be useful to other EEGLAB users with similar ' ...
               'computers as yourself.'])
    end
end
% extract result
labels = squeeze(net.getVar(net.getOutputs()).value)';
labels = reshape(mean(reshape(labels', [], 4), 2), 7, [])';