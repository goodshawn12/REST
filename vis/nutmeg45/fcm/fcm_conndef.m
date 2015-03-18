function comps=fcm_conndef(numext)
% FCM_CONNDEF defines the connections between voxel pairs to be analyzed.
%
% Usages:
%   Set the FCM configuration first with fcm_start or fcm_config_gui.
%
%   COMPS=FCM_CONNDEF;  
%   COMPS=FCM_CONNDEF(NUMEXT);
%       if you configure 'Extracerebral' channels as seed, you can optionally input
%       the number of external channels NUMEXT (e.g., number of EMG). 
% 
%   COMPS   output containing indices of voxels pairs to be analyzed.


global fuse nuts

if (~isstruct(fuse) || ~isfield(fuse,'seed'))
    error('You must set configuration first (with fcm_start or fcm_config_gui).')
end
if (~isstruct(nuts) || ~isfield(nuts,'voxels'))
    error('Open a NUTMEG session first (with nut_opensession).')
end

nv=size(nuts.voxels,1);

switch fuse.connection
    
    case 'All'
        
        switch fuse.seed
            
            case 'All'
                if isfield(fuse,'roi') && fuse.roi>0
                    nv=fuse.roi;
                end
                ncomps=(nv*(nv-1))/2;
                comps=zeros(ncomps,2,'uint32');
                currstart=1;
                for cc=1:nv
                    currlength=nv-cc;
                    currend=currstart+currlength-1;
                    comps(currstart:currend,:)=[[cc+1:nv]' cc*ones(currlength,1)];
                    currstart=currend+1;
                end
                   
            case 'Selected'    % all connections of selected voxels
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx for this FCM configuration.'), end
                nuts.selvox.ipsi=nuts.selvox.ipsi(:)';  % make sure it's a row vector
                nonsel=setdiff(1:nv,nuts.selvox.ipsi);
                ntv=length(nuts.selvox.ipsi);
                nns=length(nonsel);
                ncomps=ntv*nns + (ntv*(ntv-1))/2;
                comps=zeros(ncomps,2,'uint32');
                
                currstart=1;
                for cc=nuts.selvox.ipsi
                    curr2 = sort([nonsel(nonsel<cc) cc+1:nv])';
                    currlength=length(curr2);
                    currend=currstart+currlength-1;
                    comps(currstart:currend,:)=[cc*ones(currlength,1) curr2];
                    currstart=currend+1;
                end
                
            case 'Selected+Contralateral'
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx for this FCM configuration.'), end
                nuts.selvox.ipsi=nuts.selvox.ipsi(:)';  % make sure it's a row vector                
                if ~isfield(nuts.selvox,'contra'), nuts=fcm_contravoxidx(nuts); end
                nuts.selvox.contra=nuts.selvox.contra(:)';
                selvox=unique(cat(2,nuts.selvox.ipsi,nuts.selvox.contra));
                nonsel=setdiff(1:nv,selvox);
                ntv=length(selvox);
                nns=length(nonsel);
                ncomps=ntv*nns + (ntv*(ntv-1))/2;
                comps=zeros(ncomps,2,'uint32');

                currstart=1;
                for cc=selvox
                    curr2 = sort([nonsel(nonsel<cc) cc+1:nv])';
                    currlength=length(curr2);
                    currend=currstart+currlength-1;
                    comps(currstart:currend,:)=[cc*ones(currlength,1) curr2];
                    currstart=currend+1;
                end
                
            case 'Contralateral'
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx for this FCM configuration.'), end
                if ~isfield(nuts.selvox,'contra'), nuts=fcm_contravoxidx(nuts); end        
                nuts.selvox.contra=nuts.selvox.contra(:)';
                nonsel=setdiff(1:nv,nuts.selvox.contra);
                ntv=length(nuts.selvox.contra);
                nns=length(nonsel);
                ncomps=ntv*nns + (ntv*(ntv-1))/2;
                comps=zeros(ncomps,2,'uint32');

                currstart=1;
                for cc=nuts.selvox.contra
                    curr2 = [nonsel(nonsel<cc) cc+1:nv]';
                    currlength=length(curr2);
                    currend=currstart+currlength-1;
                    comps(currstart:currend,:)=[cc*ones(currlength,1) curr2];
                    currstart=currend+1;
                end
                
            case 'Extracerebral'
                if nargin<1, numext=1; end
                comps = uint32([]);
                for k=1:numext
                    comps = cat(1,comps,[[1:nv]' k*ones(nv,1)]);
                end
                
            otherwise
                error('This FCM configuration is not known or not programmed.')
        end
        fprintf('FCM_CONNDEF defined %d connections between voxel/ROI pairs.\n',size(comps,1))
    
    case 'Grid'   % connection of voxels with test grid
        
        switch fuse.seed
            case 'All'
                lim=[1:size(nuts.voxels,1)];
            case 'Selected'
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx first.'), end
                lim=nuts.selvox.ipsi';
            case 'Contralateral'
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx first.'), end
                lim=nuts.selvox.contra';
            case 'Selected+Contralateral'
                if ~isfield(nuts,'selvox'), error('You must run fcm_selvoxidx first.'), end
                lim = unique(cat(1,nuts.selvox.ipsi,nuts.selvox.contra))';
            otherwise
                error('Unknown seed voxel configuration.')
        end     
        grid = fcm_gridvoxidx(nuts,fuse.gridspacing);
        
        nonselgrid  = setdiff(grid,lim);
        selgrid     = intersect(lim,grid);
        selnongrid  = setdiff(lim,selgrid);
        num_sel         = length(lim);
        num_nonselgrid  = length(nonselgrid);
        num_selgrid     = length(selgrid);
        num_selnongrid  = num_sel-num_selgrid;
        num_grid        = length(grid);

        %ncomps=num_sel*num_nonselgrid + num_selgrid*num_selnongrid + (num_selgrid*(num_selgrid-1))/2;
        ncomps=num_sel*num_nonselgrid + num_selgrid*num_selnongrid + num_selgrid*(num_selgrid-1);
        comps=zeros(ncomps,2,'uint32');

        currstart=1;
        for cc=selgrid
            currend=currstart+num_grid-2;
            comps(currstart:currend,:)=[cc*ones(num_grid-1,1) setdiff(grid,cc)'];
            currstart=currend+1;
        end
        for cc=selnongrid
            currend=currstart+num_grid-1;
            comps(currstart:currend,:)=[cc*ones(num_grid,1) grid'];
            currstart=currend+1;
        end
        comps=sortrows(comps);
        
        fprintf('FCM_CONNDEF defined %d connections between voxel pairs, using a grid of %d voxels.\n',ncomps,num_grid)
            
    otherwise
        error('Unknown connection configuration.')
    
end



