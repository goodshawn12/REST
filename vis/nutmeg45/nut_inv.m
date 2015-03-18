function InvA = nut_inv(A,invtype,condthresh, param, mineigUnave)
% InvA = nut_inv(A,invtype,condthresh, param, mineigUnave)
% invtype:      choose from: tikhonov, pinv, bayesian, vanilla, or none
% condthresh:   optional, choose from: 'auto' (default) or 'always'
% param:        optional
% mineigUnave:  optional

global ndefaults 

if strcmp(invtype,'none')     % This allows easy control in functions calling nut_inv
    if(~isempty(ndefaults))
        if ndefaults.bf.viewinv
            disp(['FYI: no regularization specified (using Matlab''s ''inv'' function.)']);
        end
    end
    InvA = inv(A);
    return
end

if(nargin<3 || isempty(condthresh))
    % condthresh = 1e12; %valid to have thresh same for all types of inversion?
    condthresh = 'auto';   % changed by adrian on 16Nov2009, to avoid confusing warning message if user did not specify threshold
end

if(ischar(condthresh) && strcmp(condthresh,'always'))
    regflag = true;
else
    if( isnumeric(condthresh) && cond(A)>condthresh )
        warning('your condition number threshold kicked in, but we don''t do it this way anymore...')
    end

    % determine whether fancy regularization is necessary
    % if maximum error is greater than tenth of minimum entry
    warning('off','MATLAB:nearlySingularMatrix')
    InvA = inv(A);

    inverr = norm(inv(InvA) - A,'fro');
    warning('on','MATLAB:nearlySingularMatrix')
    minabsdata = min(abs(A(:)));
    %     if(~isempty(ndefaults))
    %         if ndefaults.bf.viewinv==1
    %             inverr
    %             minabsdata
    %         end
    %     end

    regflag = (inverr > 0.1*minabsdata ); % activate regularization in this condition
end

if(regflag)
    disp(['FYI: regularizing with ' invtype ' method.']);
    switch(invtype)
        case 'tikhonov' % param means multiplier of maximum
            if( ~exist('param','var') || isempty(param) )
                param = 1e-6;
            end
            mu = param*norm(A); % equal to param*max(eig(A))
            InvA = inv(A+mu*eye(size(A)));
        case 'tikpercent' % here, param means multiplier of diff between max and min
            if( ~exist('param','var') || isempty(param) )
                param = 1e-6;
            end
            mu = param*(max(eig(A))-min(eig(A)));
            InvA = inv(A+mu*eye(size(A)));
        case 'svdtrunc'
            warning('got here')
            [InvA]=nut_eiginv(A,1:62);
        case 'pinv'
            InvA=pinv(A);
        case 'bayesian'
            error('bayesian not ready yet.');
            gamma = var(bolts.meg,1);
            InvA = inv(A+(diag(gamma)/(size(nuts.meg.data,1)+1)));
        case 'minEigUnAve' %thanks to Matt Brookes for this version
            InvA = inv(A+(mineigUnave-min(eig(A)))*eye(size(A)));
        case 'vanilla'
            InvA = inv(A);
    end
else
    switch(invtype)
        case 'minEigUnAve' %do this way regardless
            InvA = inv(A+(mineigUnave-min(eig(A)))*eye(size(A)));
        otherwise
            if(~isempty(ndefaults))
                if ndefaults.bf.viewinv
                    disp(['FYI: no regularization needed (using Matlab''s ''inv'' function.)']);
                end
            end
            % InvA = inv(A);
            % We do not need to recalculate InvA here, because it is
            % already in workspace!
    end
end

