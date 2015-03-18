function out=fcm_fisherZtf(C,method,direc)
% FCt = fcm_fisherZtf(FC,method,direction)

if nargin<2, method='ccohere'; end
if nargin<3, direc=1; end
%if max(C(:))>1, method='Z'; end

if direc>0
    % Z-transform 
    switch method
        case {'icohere' 'glcohere'}
            out = abs(atanh(C));
        case 'ampcorr'
            out = atanh(C);  % do not need absolute values for ampcorr, as negative correlations are not merely due to ambiguity in dipole orientations
        case {'ccohere' 'nccohere'}
            C = (C./abs(C)) .* atanh(abs(C));
            % Get absolute values
            % This is necessary because the beamformer is arbitrary with regards to polarity of the virtual channel.
            try
                out = complex( abs(real(C)) , abs(imag(C)) );      
            catch ME
                if strcmp(ME.identifier,'MATLAB:nomem');
                    %out = zeros(size(C));
                    for k=1:size(C,1)
                        C(k,:,:,:) = complex( abs(real(C(k,:,:,:))) , abs(imag(C(k,:,:,:))) );
                    end
                    out=C;
                else
                    rethrow(ME)
                end
            end
        case 'pli'
            out = abs(C);
        otherwise
            out = C;
    end
elseif direc<0
    % Inverse
    switch method
        case {'ampcorr' 'icohere' 'glcohere'}
            out = tanh(C);
        case {'ccohere' 'nccohere'}
            out = (C./abs(C)) .* tanh(abs(C));
        otherwise
            out = C;
    end
end