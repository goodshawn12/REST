function elp = elec_emse2matlab(filename)

% elec_emse2matlab - Read an EMSE probe file (*.elp)
% 
% Usage: elp = elec_emse2matlab(filename)
% 
% This script extracts x,y,z values from an EMSE 
% probe (*.elp) file
% 
% EMSE *.elp files are in meters.  
% 
% Adapted for Locator v5 from Darren Webber's function under GNU GPL.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

[fid,msg] = fopen(file,'r');
if ~isempty(msg), error(msg); end

% ver = '$Revision: 1.10 $';
% fprintf('\nELEC_EMSE2MATLAB [v %s]\n',ver(11:15));
% fprintf('...reading .elp data.\n');


elp = read_elp(fid);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elp] = read_elp(fid)
    
    elp = [];
    
    % Probe files contain position information for electrode locations 
    % and/or gradiometer locations. The file consists of a prolog, a 
    % header, and a list of one or more sensor fields.
    
    % Any line beginning with '//' is a comment line, which is ignored
    
    % Read the prolog
%     tmp = fscanf(fid,'%d',2);
%     elp.version = tmp(1);
%     elp.filetype = tmp(2); % type 2 is a probe file, extension .elp
%     fgetl(fid); fgetl(fid); tmp=fgetl(fid);
%     tmp = fscanf(fid,'%d',1);
%     elp.minor_rev = tmp; % usually 1
%     
    % Read the header
    % The header consists of one optional entry and 2 entries in 
    % mandatory sequence and one optional entry:
    % Name [optional] > %N %s replace %s with name string (8 or fewer characters)
    % Type Code       > %x    replace %x with 1 (all electric), 2 (all magnetic) or 4 (mixed).
    % #Channels       > %d    number of points per channel per epoch
    
    % Sensor state (which appears in the 'type code' field) may 
    % be obtained by logically OR-ing suitable combinations from 
    % Table A-3. Note that not all combinations are physically valid.
    %
    % type/state        type code 
    % magnetic          200
    % electric          400
    % off               800
    % reference         1000
    % optical           4000
    % trigger           8000
    % other             10000
    % named point       20000
    
    % Other types (such as named points, trigger, and optical) should 
    % be represented in the same pattern as electrodes, with the type 
    % code set to identify the type. Even those types (e.g. trigger) 
    % which do not have a true location, should have a nominal 
    % location, (e.g. 0 0 0).
    
    while 1,
        tmp = fgetl(fid); % This should be: //TypeCode	nsensors
        if strmatch('//Probe type',tmp),
            tmp = fgetl(fid);
            tmp = sscanf(tmp,'%d',2);
            elp.sensorType = tmp(1);
            elp.sensorN = tmp(2);
            break;
        end
    end
    
    % Fiducial points may be included optionally. They are required 
    % for MRI registration. If they are included, they must be in 
    % the obligatory order : nasion, left preauricular point, 
    % right preauricular point. Table A-2 defines the format for 
    % representing fiduciary points.
    
    n = 0;
    while n <= 2,
        n = n + 1;
        tmp = fgetl(fid); % This should be: //Fiducials:  Nasion  Left  Right
        if strmatch('//Position of fiducials',tmp),
            tmp = fgetl(fid);
            tmp = sscanf(tmp,'%2c %f %f %f');
            elp.nasion = [tmp(3) tmp(4) tmp(5)];
            tmp = fgetl(fid);
            tmp = sscanf(tmp,'%2c %f %f %f');
            elp.lpa    = [tmp(3) tmp(4) tmp(5)];
            tmp = fgetl(fid);
            tmp = sscanf(tmp,'%2c %f %f %f');
            elp.rpa    = [tmp(3) tmp(4) tmp(5)];
            break;
        end
    end
    
    elp.x = zeros(elp.sensorN - 1,1);
    elp.y = zeros(elp.sensorN - 1,1);
    elp.z = zeros(elp.sensorN - 1,1);
    elp.ref = [];
    elp.origin = [];
    
    n = 1;
    while n <= elp.sensorN,
        
        tmp = fgetl(fid);
        if ~ischar(tmp),
            break;
        elseif strmatch('//',tmp);
            % Ignore the comment lines, get the next one
            tmp = fgetl(fid);
        end
        
        % Each electrode is represented by an electric sensor, 
        % and consists of 5 fields, of which 1 (the name) is 
        % optional. The electric sensor field data is shown 
        % in Table A-6.
        % Name              Format      Description 
        % Type Code         %S          %x replace %x with 400 (electrode) or 1c00 if reference channel
        % Name [optional]   %N          %s replace %s with name string (8 or fewer characters)
        % Position          %g %g %g    electrode location with respect to head frame (Cartesian, meters)
        % Orientation       %g %g %g    not used, replace with 0 0 1
        
        if strmatch('%S',tmp),
            
            if findstr('c00',tmp),
                ref = 1; % A reference sensor
            else
                ref = 0;
                %tmp = sscanf(tmp,'%2c %d');
                elp.type{n,1} = tmp(4:end);
            end
            
            tmp = fgetl(fid);
            if strmatch('//',tmp);
                % Ignore the comment lines, get the next one
                tmp = fgetl(fid);
            end
            
            tmp = deblank(tmp);
            if strmatch('%N',tmp),
                
                % Read the name of the sensor
                tmp = strrep(tmp,'%N','');
                tmp = fliplr(deblank(fliplr(tmp)));
                if ~ref, elp.label{n,1} = tmp; end
                
                % Read the location XYZ
                tmp = fgetl(fid);
                if strmatch('//',tmp);
                    % Ignore comments, get the next line
                    tmp = fgetl(fid);
                end
                
                if strmatch('%O',tmp),
                    if isempty(elp.origin),
                        % Get the sphere origin
                        elp.origin = sscanf(tmp(3:end),'%f',3)';
                    end
                    tmp = fgetl(fid);
                    tmp = fgetl(fid);
                    % Read the xyz location
                    tmp = sscanf(tmp,'%f',3);
                    if ref,
                        elp.ref = tmp';
                    else
                        elp.x(n) = tmp(1);
                        elp.y(n) = tmp(2);
                        elp.z(n) = tmp(3);
                        n = n + 1;
                    end
                    % Skip the next line (empty)
                    tmp = fgetl(fid);
                else
                    tmp = sscanf(tmp,'%f',3);
                    if ref,
                        elp.ref = tmp';
                    else
                        elp.x(n) = tmp(1);
                        elp.y(n) = tmp(2);
                        elp.z(n) = tmp(3);
                        n = n + 1;
                    end
                end
            end
        end
    end
    
    fclose(fid);
    
return



