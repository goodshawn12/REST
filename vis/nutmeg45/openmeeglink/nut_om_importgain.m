function [Lp,voxels] = nut_om_importgain

global nuts

ctfloadver = 'new'

omgainfiles = dir([subj '_MEGgain*.bin']);
omvoxfiles = dir([subj(1:(end-1)) '_voxels*.bin']);
if(isempty(omvoxfiles))
    omvoxfiles = dir([subj '_voxels*.bin']);
end

g=[];
voxels=[];

% join gain/voxel files
% [openmeeg calculation may have been split for memory reasons]
for ii=1:length(omgainfiles)
    g = [g om_load_full(omgainfiles(ii).name)];
    voxels = [voxels;om_load_full(omvoxfiles(ii).name)];
end

voxels = voxels(1:3:end,1:3);



% switch(nuts.meg.grad_order) % for synthetic gradient correction
%     case {1,2,3}
%         disp('importing OM lead field with synthetic gradient correction');
%     case 0
%         disp('importing OM lead field; ignoring reference channels/correction');
% end


% construct channel mixing matrix (for gradiometer derivations and
% optionally synthetic gradient correction using reference channels)
% chanmixmtx=[];
% switch(ctfloadver)
%     case 'new'
%         chanmixmtx = nuts.meg.chanmixMtx;
%     case 'old'
%         for ii=1:275
%             if(nuts.meg.grad_order) % for synthetic gradient correction
%                 chanmixmtx(ii,1:9)=-nuts.meg.Gcoef(ii,1:9);
%                 for jj=1:20
%                     chanmixmtx(ii,2*jj+[8:9]) = nuts.meg.Gcoef(ii,jj+9) * [-1 -1];
%                 end
%             end
%             chanmixmtx(ii,2*ii+[48:49])=[1 1];
%         end
% end


if(exist('ieegflag','var'))
    switch(ieegflag)
        case 1
            nuts.meg.chanmixMtx{1}=eye(size(nuts.meg.data,2));
            nuts.meg.chanmixMtx{1}(:,end+1)=-1;
        case 2
            % ignore reference electrode
            nuts.meg.chanmixMtx{1}=eye(size(nuts.meg.data,2));
            nuts.meg.chanmixMtx{1}(:,end+1)=0;
    end
end

Lp = nuts.meg.chanmixMtx{1}*g;
Lp = reshape(Lp,size(Lp,1),3,size(Lp,2)/3); % Mchannels x 3 orientations x Nvoxels

