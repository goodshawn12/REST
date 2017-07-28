% Make new eyeCatch library

% find ori_map on eyeCatch webiste, filename:
% 'eyeChannelWeightNormalizedpart1.mat' and 'eyeChannelWeightNormalizedpart1.mat'

function new_map = interpo_lib(ori_map, downsize)

[sampleNumber, oldsize] = size(ori_map);
if downsize > sqrt(oldsize)
    error('Downsize map should be smaller than original map.');
end

oldsize = sqrt(oldsize);

% reshape old_map
ori_map = reshape(ori_map, [sampleNumber, oldsize, oldsize]);
new_map = zeros(sampleNumber, downsize^2);

%% resize
for i = 1:sampleNumber
    new_map(i,:) = reshape(imresize(squeeze(ori_map(i,:,:)),[32 32],'nearest'), [1, downsize^2]);
end

% save to mat
filename = 'libEyeCatch';

% for oversize library
% k = 1;
% filename = ['new_lib',int2str(k)];
% while exist([filename,'.mat'],'file')
%     k = k+1;
%     filename = ['new_lib.mat', int2str(k)];
% end

save(filename, 'new_map');

end




