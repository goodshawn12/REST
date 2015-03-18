function origin_mm = nut_read_mriorigin(fn,convention)
% NUT_READ_MRIORIGIN  finds origin of structural MRI in mm
%
%   origin_mm = nut_read_mriorigin('mrifile.hdr',convention)
%
% mrifile.hdr       Name of analyze header file
% convention        [optional] 'SPM' (default) or 'SMAC'

if nargin<2
    convention='SPM';
end

switch upper(convention)
    case 'SPM'
        
        % adapted from SPM
        hdr=spm_read_hdr(fn);
        if any(hdr.hist.origin(1:3)),
            origin_pix = hdr.hist.origin(1:3);
        else
            origin_pix = (hdr.dime.dim(2:4)+1)/2;
        end
        origin_mm = origin_pix .* hdr.dime.pixdim(2:4);
        
    case 'SMAC'
        
        if isempty(which('mriCalcHistogram'))
            error('NUTMEG:noSMACsrc','Could not find SMAC m-files in your MATLAB path.')
        end
        
        % adapted from SMAC
        [mri,hdr] = ioReadMRI(fn);
        mriProperties.mriHistogram = mriCalcHistogram(mri, max(mri(:)), '_pro');
        mriProperties.thresholdValue = cast(max(1, mriProperties.mriHistogram.bgLevel), 'double');
        mriProperties.boundingBox = mriGetBoundingBox(mri, mriProperties.thresholdValue);
        origin_pix = round(sum(mriProperties.boundingBox,1)./2);
        voxdim = hdr.image_dimension.pixdim(2:4);
        if ~isequal(voxdim, [1 1 1])
            warning('NUTMEG:badMRIvoxelsize','Your MRI does not seem to have a sampling rate of 1x1x1 mm. Make sure you have the right voxelsize.')
            origin_mm = origin_pix .* voxdim;
        else
            origin_mm = origin_pix;
        end
end
