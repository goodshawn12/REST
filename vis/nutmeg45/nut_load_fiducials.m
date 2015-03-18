function nut_load_fiducials(fullpath_lsc,V)
% function nut_load_fiducials(fullpath_lsc,V)
%
% fullpath_lsc: full path/file name of header file with fiducial
%               information (e.g. *.hdm file for CTF)
% V (optional): if calling from command line outside of Nutmeg, sets which
%               MRI to which the fiducials are alligned. V is SPM result
%               from calling V=spm_vol(mripath); However, if Nutmeg open
%               then this will be in global st (st.vols{1})

global coreg st
if exist('V','var')
    st.vols{1}=V;
end
if isempty(st)
    error('you need some mri info either input V or from global st')
end
if exist('fullpath_lsc','var')
    [pathname_lsc,filename_lsc,filetype]=fileparts(fullpath_lsc);
    filename_lsc=[filename_lsc filetype];
else
    [filename_lsc, pathname_lsc]=uigetfile('*.txt;*.hdm','Please Select the stored fiducials (CTF .hdm file supported)');
    if isequal(filename_lsc,0)|isequal(pathname_lsc,0)
        return;
    end
    filetype = filename_lsc(end-3:end);
end
switch filetype
    case '.txt'
        
        [x y z] =textread([pathname_lsc filename_lsc],'%n%n%n%*[^\n]',3);
        coreg.fiducials_mri_mm=[x y z];
        if regexpi(st.vols{1}.fname,'blank.img')
            %nuts.coreg=nut_CoregistrationTool(nuts.coreg);
            %FIXME: create button to load just lsc in absence of MRI
            disp('you need an MRI before loading fiducials.  you may load LSC without MRI');
        else
            coreg = nut_coreg_fiducials(false,coreg);
        end


%         hx_obj = findobj('Tag','nut_left_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(1),y(1),z(1)));
%         % Switch the "Show left PA" button ON
%         %set(findobj('Tag','nut_showleft'),'Enable','on')
%         hx_obj = findobj('Tag','nut_right_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(2),y(2),z(2)));
%         % Switch the "Show right PA" button ON
%         %set(findobj('Tag','nut_showright'),'Enable','on')
%         hx_obj = findobj('Tag','nut_nose_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(3),y(3),z(3)));
%         % Switch the "Show Nasion" button ON
%         %set(findobj('Tag','nut_shownose'),'Enable','on')
    case '.hdm'    %%% get lsc from CTF's head model file

        lscfid=fopen(fullfile(pathname_lsc,filename_lsc));

        % get the local sphere center(s) -- stored in MEG coords
        while(1)
            temp = fgetl(lscfid);
            % if EOF reached and multisphere not found, go back and get
            % single lsc
            if(~ischar(temp))
                fseek(lscfid,lscseek,'bof');
                fgetl(lscfid); % skip next line with '{'
                fscanf(lscfid,'%s',1); % should be "ORIGIN_X:"
                x = fscanf(lscfid,'%f',1);
                fscanf(lscfid,'%s',1); % should be "ORIGIN_Y:"
                y = fscanf(lscfid,'%f',1);
                fscanf(lscfid,'%s',1); % should be "ORIGIN_Z:"
                z = fscanf(lscfid,'%f',1);
                fscanf(lscfid,'%s',1); % should be "RADIUS:"
                r_head = fscanf(lscfid,'%f',1);   
                xyzr = [x y z r_head] *10; % given in cm, we need mm
                
                %hx_obj = findobj('Tag','nut_lsc_text');
                %set the static text box with the above value
                %set(hx_obj,'String',sprintf('MEG: %.1f %.1f %.1f',x,y,z));
                break; % bust out of while loop
            end

            %  record single lsc position in case multisphere not present
            if(strcmp(temp,'MEG_Sphere'))
                lscseek = ftell(lscfid);
            end

            % multiple sphere data...
            if(strcmp(temp,'MultiSphere_Data'))
                % skip next 7 lines
                for ii=1:7
                    fgetl(lscfid);
                end
                testname = fscanf(lscfid,'%s',1);
                if strcmp(testname,'HEADPOS:')
                    % skip next 3 lines
                    disp('CTF_HEAD_MODEL_FILE_VERSION_6.0 detected');
                    for ii=1:3
                        fgetl(lscfid);
                    end
                    clear testname
                else
                    channelname=testname;
                    clear testname
                end

                channelnum=0;
                while(1)
                    if ~exist('channelname') % it will exist if testname above is a channel
                        channelname = fscanf(lscfid,'%s',1);
                    end
                    if(strcmp(channelname,'}'))
                        break; % end of channel list, bust out of while loop
                    end
                    channelnum = channelnum + 1;
                    lsc_sensor_labels{channelnum} = channelname(1,1:(end-1));
                    xyzr(channelnum,:) = fscanf(lscfid,'%f',4);
                    clear channelname
                    
                end
                xyzr=xyzr*10;  % given in cm, we need mm
                r_head = xyzr(:,4);
                
%                 hx_obj = findobj('Tag','nut_lsc_text');
%                 %set the static text box with the above value
%                 set(hx_obj,'String','multiple spheres');
                
                break % bust out of while loop
            end
        end

        % insert lsc info into nuts.meg
        global nuts
        %nuts.meg.lsc = 10*(xyzr(:,1:3));
        for ii=1:3,  %% lsc can't be integer since will later difference with voxel coords, and will produce divide by zero error
            if(round(xyzr(:,ii))==xyzr(:,ii))
                xyzr(:,ii)=xyzr(:,ii)+0.01;
                disp('adding 0.01 to lsc coord to avoid divide by zero issue later on')
            end
        end
        nuts.meg.lsc = (xyzr(:,1:3));
        if(exist('lsc_sensor_labels','var'))
            nuts.meg.lsc_sensor_labels = lsc_sensor_labels;
        elseif(isfield(nuts.meg,'lsc_sensor_labels'));
            nuts.meg = rmfield(nuts.meg,'lsc_sensor_labels');
        end
        clear nuts

        % now jump down to fiducial block
        fseek(lscfid,0,'bof');
        tmpstr=fgetl(lscfid);
        while(~any(strcmp(tmpstr,{'Fid_Points';'MRI_Fid_Points'})))
            if(tmpstr==-1)
                disp('this is weird, you have sphere center(s) in your default.hdm but not fiducials?  okay, then using the ones already in your nuts.coreg')
                return;
            end
            tmpstr=fgetl(lscfid);
        end
        fgetl(lscfid); % skip next line with "{"
        fgetl(lscfid); % skip next line with "// Fid. Pts"
        fscanf(lscfid,'%s',1); % should be "NASION:"
        n = fscanf(lscfid,'%f',3); % nasion
        fscanf(lscfid,'%s',1); % should be "LEFT_EAR:"
        l = fscanf(lscfid,'%f',3); % left
        fscanf(lscfid,'%s',1); % should be "RIGHT_EAR:"
        r = fscanf(lscfid,'%f',3); % right
        fclose(lscfid);
        
        lrn = [l r n]';
        
        if ~isempty(st) && ( ~isempty(st.vols{1}) && max(size(regexpi(st.vols{1}.fname,'/blank.img$'))) )
            %nuts.coreg=nut_CoregistrationTool(nuts.coreg);
            %FIXME: create button to load just lsc in absence of MRI
            disp('you need an MRI before loading fiducials.  you may load LSC without MRI');
        elseif (isempty(coreg)) %this will assumedly happen when being called from nut_results_viewer
            
        else
            DIM = st.vols{1}.dim(1:3);
            VOX = nut_pixdim(st.vols{1});
            
            % image orientation affects fiducial coord conversion, so check for it here:
            % value of 1 = neurological (L is L); 2 = radiological (L is R)
            %coreg.orientation = get(findobj('Tag','nut_orientation_menu'),'Value');

            % convert from CTF's isotropic coords to standard Analyze coords
            % assumes 256x256x256 CTF volume and 256x256xXXX Analyze volume...
            
            if(DIM(1) ~= 256 || DIM(2) ~= 256)
                warndlg('Your MRI is a non-standard size... I don''t know how to properly convert the HDM file. Do not trust these fiducials!')
            end
    
            % if radiological, flip x-coords
            if ~isfield(coreg,'orientation')
                coreg.orientation=1;
            end
            if(~isempty(coreg.orientation) && coreg.orientation == 2)
                ctf2spm_tfm = [1 0 0 0; 0 1 0 257; 0 0 -VOX(1)/VOX(3) [DIM(3)/2 + 128*VOX(1)/VOX(3) + 1]; 0 0 0 1];
            else
                ctf2spm_tfm = [1 0 0 0; 0 -1 0 257; 0 0 -VOX(1)/VOX(3) [DIM(3)/2 + 128*VOX(1)/VOX(3) + 1]; 0 0 0 1];
            end
            xyz_vx = nut_coordtfm(lrn, ctf2spm_tfm);
            coreg.fiducials_mri_mm=nut_voxels2mm(xyz_vx);
%%% OLD WAYS, probably more wrong.
%             x_vx=lrn(:,1);
%             if(~isempty(coreg.orientation) && coreg.orientation == 2)  % if radiological
%                 x_vx = 257-x_vx;  % flip x-coords as well
%             end
%             y_vx=257-lrn(:,2);    % y-coords always need to be flipped
%             z_vx=(0.5*(256+DIM(3)*VOX(3)/VOX(1))-lrn(:,3))*VOX(1)/VOX(3) +1 + VOX(3)/VOX(1);
%             %z_vx=(0.5*(256+DIM(3)*VOX(3)/VOX(1))-lrn(:,3))*VOX(1)/VOX(3);  %older way
%             [x y z]=nut_voxels2mm(x_vx,y_vx,z_vx);
%             coreg.fiducials_mri_mm=[x y z];

            coreg = nut_coreg_fiducials(false,coreg)

        end

        
%         hx_obj = findobj('Tag','nut_left_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(1),y(1),z(1)));
% 
%         hx_obj = findobj('Tag','nut_right_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(2),y(2),z(2)));
% 
%         hx_obj = findobj('Tag','nut_nose_text');
%         %set the static text box with the above value
%         set(hx_obj,'String',sprintf('MRI: %.1f %.1f %.1f',x(3),y(3),z(3)));
        
    otherwise
        errordlg('Unsupported file type.')
end

global nuts
if ~isempty(nuts)
    if isfield(nuts,'coreg') & (isfield(nuts.coreg,'norm_mripath') & ~isempty(coreg))
        coreg.fiducials_mni_mm=nut_mri2mni(coreg.fiducials_mri_mm);
    end
else
    warning('Did you know you have no nuts?')
end
clear nuts
