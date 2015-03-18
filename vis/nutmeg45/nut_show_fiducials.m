function nut_show_fiducials
% NUT_SHOW_FIDUCIALS
%
% Show the fiducials in 3D
%
% The function uses NUTS global.

global coreg
switch get(gcbo,'Tag')
    case {'nut_show_fiducials_3d'}
        %check if the required file and variables exist
        if isempty(coreg)
         msgbox('Load the image','Volume not found','warn')
         return;
        end
        if ~isfield(coreg,'fiducials_mri_mm')
            errordlg('Dude, where''s your fiducials?');
            return;            
        end
        if   ~isfield(coreg,'mesh')
         errordlg('Dude, where''s your head?');
         return            
        end
        fiducials=coreg.fiducials_mri_mm;
        %display fiducials in 3D
        figure
        nut_show_head(coreg.mesh);
        hold on;       
        plot3(fiducials(:,1),fiducials(:,2),fiducials(:,3),'.');
        hold off
        clear fiducials
        return                
    otherwise
        errordlg('Callback not yet defined for this button',...
			'UNDEFINED CALLBACK');
        return
end



