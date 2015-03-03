% Defines the class currentSourceViewer for visualization of EEG inverse solutions on a realistic head model. 
% This class is part of MoBILAB software. 
% For more details visit:  https://code.google.com/p/mobilab/
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2012

classdef currentSourceViewer < handle
    properties
        hFigure
        hAxes
        streamObj
        hLabels
        hSensors
        hScalp
        hCortex
        hVector
        dcmHandle
        sourceMagnitud
        sourceOrientation
        scalpData
        pointer
    end
    methods
        function obj = currentSourceViewer(streamObj,J,V,figureTitle,channelLabels)
            if nargin < 3, V = [];end
            if nargin < 4, figureTitle = '';end
            if nargin < 5
                channelLabels = cell(size(V,1),1);
                for it=1:length(channelLabels), channelLabels{it} = num2str(it);end
            end
            if isempty(channelLabels)
                channelLabels = cell(size(V,1),1);
                for it=1:length(channelLabels), channelLabels{it} = num2str(it);end
            end
            obj.streamObj = streamObj; 
            load(obj.streamObj.surfaces);
            color = [0.93 0.96 1];
            
            mPath = which('mobilabApplication');
            if isempty(mPath)
                 path = fileparts(which('currentSourceViewer'));
            else path = fullfile(fileparts(mPath),'skin');
            end
            
            try labelsOn  = imread([path filesep 'labelsOn.png']);
                labelsOff = imread([path filesep 'labelsOff.png']);
                sensorsOn = imread([path filesep 'sensorsOn.png']);
                sensorsOff = imread([path filesep 'sensorsOff.png']);
                scalpOn = imread([path filesep 'scalpOn.png']);
                scalpOff = imread([path filesep 'scalpOff.png']);
                vectorOn = imread([path filesep 'vectorOn.png']);
                vectorOff = imread([path filesep 'vectorOff.png']);
                prev = imread([path filesep '32px-Gnome-media-seek-backward.svg.png']);
                next = imread([path filesep '32px-Gnome-media-seek-forward.svg.png']);
            catch ME
                ME.message = sprintf('%s\nSome icons may be missing.',ME.message);
                ME.rethrow;
            end
            if isa(streamObj,'struct'), visible = 'off';else visible = 'on';end
            obj.hFigure = figure('Menubar','figure','ToolBar','figure','renderer','opengl','Visible',visible,'Color',color,'Name',figureTitle);
            position = get(obj.hFigure,'Position');
            set(obj.hFigure,'Position',[position(1:2) 1.06*position(3:4)]);
            obj.hAxes = axes('Parent',obj.hFigure);         
            toolbarHandle = findall(obj.hFigure,'Type','uitoolbar');
            
            hcb(1) = uitoggletool(toolbarHandle,'CData',labelsOff,'Separator','on','HandleVisibility','off','TooltipString','Labels On/Off','userData',{labelsOn,labelsOff},'State','off');
            set(hcb(1),'OnCallback',@(src,event)rePaint(obj,hcb(1),'labelsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(1),'labelsOff'));
            
            hcb(2) = uitoggletool(toolbarHandle,'CData',sensorsOff,'Separator','off','HandleVisibility','off','TooltipString','Sensors On/Off','userData',{sensorsOn,sensorsOff},'State','off');
            set(hcb(2),'OnCallback',@(src,event)rePaint(obj,hcb(2),'sensorsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(2),'sensorsOff'));
            
            hcb(3) = uitoggletool(toolbarHandle,'CData',scalpOff,'Separator','off','HandleVisibility','off','TooltipString','Scalp On/Off','userData',{scalpOn,scalpOff},'State','off');
            set(hcb(3),'OnCallback',@(src,event)rePaint(obj,hcb(3),'scalpOn'),'OffCallback',@(src, event)rePaint(obj,hcb(3),'scalpOff'));
            
            hcb(4) = uitoggletool(toolbarHandle,'CData',vectorOff,'Separator','off','HandleVisibility','off','TooltipString','Vectors On/Off','userData',{vectorOn,vectorOff},'State','off');
            set(hcb(4),'OnCallback',@(src,event)rePaint(obj,hcb(4),'vectorOn'),'OffCallback',@(src, event)rePaint(obj,hcb(4),'vectorOff'));
            
            uipushtool(toolbarHandle,'CData',prev,'Separator','off','HandleVisibility','off','TooltipString','Previous','ClickedCallback',@obj.prev);
            uipushtool(toolbarHandle,'CData',next,'Separator','off','HandleVisibility','off','TooltipString','Next','ClickedCallback',@obj.next);
            set(obj.hFigure,'WindowScrollWheelFcn',@(src, event)mouseMove(obj,[], event));
            
            obj.dcmHandle = datacursormode(obj.hFigure);
            obj.dcmHandle.SnapToDataVertex = 'off';
            set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            obj.dcmHandle.Enable = 'off';
            hold(obj.hAxes,'on');
            
            obj.hSensors = scatter3(obj.hAxes,obj.streamObj.channelSpace(:,1),obj.streamObj.channelSpace(:,2),...
                obj.streamObj.channelSpace(:,3),'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
            set(obj.hSensors,'Visible','off');
            
            N = length(channelLabels);
            k = 1.1;
            obj.hLabels = zeros(N,1);
            for it=1:N, obj.hLabels(it) = text('Position',k*obj.streamObj.channelSpace(it,:),'String',channelLabels{it},'Parent',obj.hAxes);end
            set(obj.hLabels,'Visible','off');
            
            if size(J,1) == 3*size(surfData(3).vertices,1)  
                J = reshape(J,[size(J,1)/3 3 size(J,2)]);
                Jm = squeeze(sqrt(sum(J.^2,2)));
                normals = J;
            else
                Jm = J;
                normals = geometricTools.getSurfaceNormals(surfData(3).vertices,surfData(3).faces,false);    
            end
            obj.sourceMagnitud = Jm;
            obj.sourceOrientation = J;
            obj.pointer = 1;
            
            % vectors
            obj.hVector = quiver3(surfData(3).vertices(:,1),surfData(3).vertices(:,2),surfData(3).vertices(:,3),normals(:,1,1),normals(:,2,1),normals(:,3,1),2);
            set(obj.hVector,'Color','k','Visible','off');
            
            % cortex
            obj.hCortex = patch('vertices',surfData(3).vertices,'faces',surfData(3).faces,'FaceVertexCData',obj.sourceMagnitud(:,1),...
                'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
                'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes);
            camlight(0,180)
            camlight(0,0)
            mx = 0.9*max(abs(Jm(:)));
            
            % scalp
            if isempty(V)
                skinColor = [1,.75,.65];
                obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'facecolor',skinColor,...
                    'facelighting','phong','LineStyle','none','FaceAlpha',.85,'Parent',obj.hAxes,'Visible','off');
                obj.scalpData = [];
            else
                W = geometricTools.localGaussianInterpolator(obj.streamObj.channelSpace,surfData(1).vertices,6);
                obj.scalpData = W*V;
                obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'FaceVertexCData',obj.scalpData(:,1),...
                    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',0.85,'SpecularColorReflectance',0,...
                    'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes,'Visible','off');
            end
            view(obj.hAxes,[90 0]);
            colorbar
            % box on;
            hold(obj.hAxes,'off');
            axis(obj.hAxes,'equal','vis3d');
            axis(obj.hAxes,'off')
            set(obj.hAxes,'Clim',[-mx mx]);
            colormap(bipolar(512, 0.99))
            if isprop(streamObj,'name'), objName = [streamObj.name ': '];else objName = '';end
            set(obj.hFigure,'Visible',visible,'userData',obj,'Name',[objName num2str(1) '/' num2str(size(obj.sourceMagnitud,2))]);
            rotate3d
            drawnow
        end
        %%
        function rePaint(obj,hObject,opt)
            CData = get(hObject,'userData');
            if isempty(strfind(opt,'Off'))
                set(hObject,'CData',CData{2});
            else
                set(hObject,'CData',CData{1});
            end
            switch opt
                case 'labelsOn'
                    set(obj.hLabels,'Visible','on');
                case 'labelsOff'
                    set(obj.hLabels,'Visible','off');
                case 'sensorsOn'
                    set(obj.hSensors,'Visible','on');
                case 'sensorsOff'
                    set(obj.hSensors,'Visible','off');
                case 'scalpOn'
                    val = obj.scalpData(:);
                    mx = 0.9*max(abs(val));
                    set(obj.hCortex,'FaceAlpha',0.15);
                    if strcmp(get(obj.hVector,'Visible'),'on')
                        set(obj.hScalp,'Visible','on','FaceAlpha',0.65);
                    else
                        set(obj.hScalp,'Visible','on','FaceAlpha',0.85);
                    end
                    if isempty(val), return;end
                    set(get(obj.hScalp,'Parent'),'Clim',[-mx mx]);
                case 'scalpOff'
                    val = obj.sourceMagnitud(:);
                    mx = 0.9*max(abs(val));
                    set(obj.hScalp,'Visible','off');
                    set(obj.hCortex,'Visible','on','FaceAlpha',1);
                    set(get(obj.hCortex,'Parent'),'Clim',[-mx mx]);
                case 'vectorOn'
                    set(obj.hVector,'Visible','on');
                    set(obj.hCortex,'FaceAlpha',0.75);
                case 'vectorOff'
                    set(obj.hVector,'Visible','off');
                    set(obj.hCortex,'FaceAlpha',1);
            end
        end
        %%
        function output_txt = showLabel(obj,event_obj)
            persistent DT
            if strcmp(obj.dcmHandle.Enable,'off'),return;end
            if strcmp(get(obj.hVector,'Visible'),'on')
                set(obj.hVector,'Visible','off');
                set(obj.hCortex,'FaceAlpha',1);
            end
            if isempty(DT)
                load(obj.streamObj.surfaces);
                vertices = surfData(3).vertices;
                DT = DelaunayTri(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            pos = get(event_obj,'Position');
            loc = nearestNeighbor(DT, pos);
            output_txt = obj.streamObj.atlas.label{obj.streamObj.atlas.colorTable(loc)};
            drawnow
            %updateCursor(obj.dcmHandle,pos);
        end
        %%
        function prev(obj,~,~)
            obj.pointer = obj.pointer-1;
            if obj.pointer < 1, obj.pointer = 1;end
            val = obj.sourceMagnitud(:,obj.pointer);
            set(obj.hCortex,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hCortex,'Visible'),'on'), set(get(obj.hCortex,'Parent'),'Clim',[-mx mx]);end
            if isprop(obj.streamObj,'name'), objName = [obj.streamObj.name ': '];else objName = '';end
            set(obj.hFigure,'Name',[objName num2str(obj.pointer) '/' num2str(size(obj.sourceMagnitud,2))]);
            if isempty(obj.scalpData), drawnow;return;end
            val = obj.scalpData(:,obj.pointer);
            set(obj.hScalp,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hScalp,'Visible'),'on'), set(get(obj.hScalp,'Parent'),'Clim',[-mx mx]);end
            try %#ok
                set(obj.hVector,'UData',obj.sourceOrientation(:,1,obj.pointer),'VData',obj.sourceOrientation(:,2,obj.pointer),'WData',obj.sourceOrientation(:,3,obj.pointer));
            end
        end
        function next(obj,~,~)
            obj.pointer = obj.pointer+1;
            n = size(obj.sourceMagnitud,2);
            if obj.pointer > n, obj.pointer = n;end
            val = obj.sourceMagnitud(:,obj.pointer);
            set(obj.hCortex,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hCortex,'Visible'),'on'), set(get(obj.hCortex,'Parent'),'Clim',[-mx mx]);end
            if isprop(obj.streamObj,'name'), objName = [obj.streamObj.name ': '];else objName = '';end
            set(obj.hFigure,'Name',[objName num2str(obj.pointer) '/' num2str(size(obj.sourceMagnitud,2))]);
            if isempty(obj.scalpData), drawnow;return;end
            val = obj.scalpData(:,obj.pointer);
            set(obj.hScalp,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hScalp,'Visible'),'on'), set(get(obj.hScalp,'Parent'),'Clim',[-mx mx]);end
            try %#ok
                set(obj.hVector,'UData',obj.sourceOrientation(:,1,obj.pointer),'VData',obj.sourceOrientation(:,2,obj.pointer),'WData',obj.sourceOrientation(:,3,obj.pointer));
            end
            drawnow
        end
        function mouseMove(obj,~,eventObj)
            obj.pointer = obj.pointer - eventObj.VerticalScrollCount;%*eventObj.VerticalScrollAmount;
            if obj.pointer < 1, obj.pointer = 1;end
            if obj.pointer > size(obj.sourceMagnitud,2), obj.pointer = size(obj.sourceMagnitud,2);end
            val = obj.sourceMagnitud(:,obj.pointer);
            if isprop(obj.streamObj,'name'), objName = [obj.streamObj.name ': '];else objName = '';end
            set(obj.hFigure,'Name',[objName num2str(obj.pointer) '/' num2str(size(obj.sourceMagnitud,2))]);
            set(obj.hCortex,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hCortex,'Visible'),'on'), set(get(obj.hCortex,'Parent'),'Clim',[-mx mx]);end
            if isempty(obj.scalpData), drawnow;return;end
            val = obj.scalpData(:,obj.pointer);
            set(obj.hScalp,'FaceVertexCData',val);
            mx = 0.9*max(abs(val));
            %if strcmp(get(obj.hScalp,'Visible'),'on'), set(get(obj.hScalp,'Parent'),'Clim',[-mx mx]);end
        end
    end
end
