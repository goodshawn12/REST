% The static class geometricTools encapsulates the most used functions for surface processing
% used to solve EEG forward or inverse problems.
%
% Author: Alejandro Ojeda, SCCN/INC/UCSD, 2012
%
% Contributors:
% nonrigid_version23    -> D.Kroon, University of Twente, August
%                             2010 (http://www.mathworks.com/matlabcentral/fileexchange/20057)
% getSurfaceLaplacian   -> Nelson Trujillo Barreto and Pedro Antonio Valdes
%                             Hernandez, Cuban Neuroscience Center
% iso2mesh dependencies -> Qianqian Fang, http://iso2mesh.sourceforge.net/cgi-bin/index.cgi

classdef geometricTools
    methods
        function obj = geometricTools()
        end
        %%
        function disp(obj)
            disp(['Static class: ' class(obj)])
            disp('Static classes work like a namespace, they are used to create data and functions that can be accessed without creating an instance of the class.')
            methods(obj);
        end
    end
    methods(Static)
        %%
        function Xcentered = correctOrigin(X)
            [m,n] = size(X);
            K = ones(m,n);
            B = pinv(K'*K)*K'*X;
            X0 = K*B;
            Aff = eye(4);
            Aff([1 2],4) = X0(1,[1 2])';
            Xcentered = Aff\[X ones(m,1)]';
            Xcentered = Xcentered(1:3,:)';
        end
        %%
        function [Aff,Sn, scale] = affineMapping(S,T)
            % S: source space
            % T: target space
            % S = [sx1 sy1 sz1; sx2 sy2 sz2; ... sxk syn szk]
            % T = [tx1 ty1 tz1; tx2 ty2 tz2; ... txk tyn tzk]
            %
            % d(Aff) = min frobenius(T -S* Aff')
            % Sn = S*Aff';
            
            [~,~,transform] = procrustes(T,S);
            scale = transform.b;
            Aff = [[transform.b*transform.T;transform.c(1,:)] [0;0;0;1]]';
            Sn = geometricTools.applyAffineMapping(S,Aff);
        end
        %%
        function T = applyAffineMapping(S,M)
            % [T 1] = [S 1]*M';
            T = [S ones(size(S,1),1)]*M';
            T(:,4) = [];
        end
        %%
        function [def,spacing,offset,SgridWarped] = bSplineMapping(S,T,Sgrid,options)
            if nargin < 4,
                options.Verbose = false;
                options.MaxRef = 5;
            end
            mn = min(Sgrid);
            Smn = bsxfun(@minus,S,mn);
            dim = max(Sgrid) - mn;
            Tmn = bsxfun(@minus,T,mn);
            [def,spacing,SgridWarped] = point_registration(dim,Smn,Tmn,options);
            offset = mn;
            SgridWarped = bsxfun(@plus,SgridWarped,offset);
        end
        %%
        function SgridWarped = applyBSplineMapping(def,spacing,offset,Sgrid)
            Smn = bsxfun(@minus,Sgrid,offset);
            SmnWarped = bspline_trans_points_double(def,spacing,Smn);
            SgridWarped = bsxfun(@plus,SmnWarped,offset);
        end
        %%
        function [neighbors,D,loc] = nearestNeighbor(vertices,T)
            if exist('DelaunayTri','file')
                 dt = DelaunayTri(vertices(:,1),vertices(:,2),vertices(:,3)); %#ok
            else dt = delaunayTriangulation(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            try [loc,D] = nearestNeighbor(dt, T);
            catch
                loc = nearestNeighbor(dt, T);
                D = sqrt(sum((vertices(loc,:)-T).^2,2));
            end
            neighbors = vertices(loc,:);
        end
        %%
        function Yi = ridgeInterpolation(vertices,faces,elec,Y)
            L = geometricTools.getSurfaceLaplacian(vertices,faces);
            K = geometricTools.localGaussianInterpolator(vertices,elec,1);
            K = full(K);
            Yi = ridgeGCV(Y,K,L,100,0);
        end
        %%
        function J = simulateGaussianSource(X,X0,h)
            if nargin < 3, h = 0.1;end
            J = geometricTools.localGaussianInterpolator(X,X0,h)';
        end
        function W = localGaussianInterpolator(X,Xi,h,normalize)
            if nargin < 3, h = 0.1;end
            if nargin < 4, normalize = false;end
            N = size(Xi,1);
            M = size(X,1);
            W = zeros(N,M);
            for it=1:N
                d = sum(bsxfun(@minus,X,Xi(it,:)).^2,2);
                W(it,:) = exp(-d/(2*h^2));
            end
            if normalize, W = bsxfun(@rdivide,W,sum(W,2)+eps);end
        end
        %%
        function D = isInConvexHull(X,Xi)
            
            X = geometricTools.correctOrigin(X);
            X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
            
            [x,y,z] = sphere(64);
            
            figure;
            surf(x,y,z,'FaceColor','g','FaceAlpha',0.7,'EdgeColor','none');
            set(gca,'Projection','perspective','DataAspectRatio',[1 1 1]); hold on;axis tight;camlight
            
            plot3(X(:,1),X(:,2),X(:,3),'.')
            
            X = geometricTools.correctOrigin(X);
            X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
            
            
            N = size(Xi,1);
            M = size(X,1);
            D = zeros(M,N);
            for it=1:N
                d = bsxfun(@minus,X,Xi(it,:));
                D(:,it) = sqrt(sum(( d ).^2,2));
            end
        end
        %%
        function Yi = localGaussianScatterInterpolator(X,Y,Xi)
            n = size(Y,2);
            if n==1
                F = TriScatteredInterp(X(:,1),X(:,2),X(:,3),Y,'nearest');
                Yi = F(Xi(:,1),Xi(:,2),Xi(:,3));
            else
                Yi = zeros(size(Xi,1),n);
                for it=1:n
                    F = TriScatteredInterp(X(:,1),X(:,2),X(:,3),Y(:,it),'nearest');
                    Yi(:,it) = F(Xi(:,1),Xi(:,2),Xi(:,3));
                end
            end
            W = geometricTools.localGaussianInterpolator(Xi,Xi,16);
            Yi = W*Yi;
        end
        %%
        function [rVertices,rFaces] = resampleSurface(vertices,faces,decimationPercent)
            if nargin < 2, error('Not enough input arguments.');end
            if nargin < 3, decimationPercent = 0.1;end
            if isempty(which('meshresample')), error('This function uses Iso2Mesh toolbox, you can download it for free fom: http://iso2mesh.sourceforge.net');end
            [rVertices,rFaces]=meshresample(vertices,faces,decimationPercent);
        end
        %%
        function sVertices = smoothSurface(vertices,faces,lambda,method)
            if nargin < 2, error('Not enough input arguments.');end
            if nargin < 3, lambda = 0.2;end
            if nargin < 4, method = 'lowpass';end
            
            maxIter = 20;
            N = size(vertices,1);
            
            if isempty(which('meshresample'))
                warning('MoBILAB:noIso2Mesh','This function uses Iso2Mesh toolbox if is installed, you can download it for free fom: http://iso2mesh.sourceforge.net');
                sVertices = vertices;
                for it=1:N
                    ind = any(faces==it,2);
                    indices = faces(ind,:);
                    indices = indices(:);
                    indices(indices==it) = [];
                    W = geometricTools.localGaussianInterpolator(vertices(indices,:),vertices(it,:),10);
                    sVertices(it,:) = sum((1./W)*vertices(indices,:),1)./sum(1./W);
                end
                return;
            end
            conn = neighborelem(faces,size(vertices,1));
            for it=1:N
                tmp = faces(conn{it},:);
                conn{it} = unique(tmp(:)');
            end
            sVertices = smoothsurf(vertices,[],conn,maxIter,lambda,method);
        end
        %%
        function [rVertices,rFaces] = refineSurface(vertices,faces,decimationRate,maxIter)
            if nargin < 3, decimationRate = 0.5;end
            if nargin < 4, maxIter = 3;end
            if isempty(which('meshresample')), error('This function uses Iso2Mesh toolbox, you can download it for free fom: http://iso2mesh.sourceforge.net');end
            
            tmpVertices = vertices;
            tmpFaces = faces;
            for it=1:maxIter
                [tmpVertices,tmpFaces] = geometricTools.resampleSurface(tmpVertices,tmpFaces,decimationRate);
                tmpVertices = geometricTools.smoothSurface(tmpVertices,tmpFaces);
                disp(it)
            end
            rVertices = tmpVertices;
            rFaces = tmpFaces;
        end
        %%
        function verticesExt = repareIntersectedSurface(surfInt,surfOut,dmax)
            if nargin < 3, dmax = 8;end
            verticesInt = surfInt.vertices;
            verticesExt = surfOut.vertices;
            [nVerticesInt,d] = geometricTools.nearestNeighbor(verticesExt,verticesInt);
            I = d < dmax;
            while any(I)
                I2 = ismember(verticesExt,nVerticesInt(I,:),'rows');
                verticesExt(I2,:) = 1.005*verticesExt(I2,:);
                [nVerticesInt,d] = geometricTools.nearestNeighbor(verticesExt,verticesInt);
                I = d < dmax;
            end
            if any(verticesExt(:) ~= surfOut.vertices(:))
                verticesExt = geometricTools.smoothSurface(verticesExt,surfOut.faces);
            end
        end
        %%
        function [normals,faces] = getSurfaceNormals(vertices,faces,normalsIn)
            if nargin < 3, normalsIn = true;end
            h = figure('visible','off');
            h2 = patch('vertices',vertices,'faces',fliplr(faces));
            normals = get(h2,'vertexnormals');close(h);
            if isempty(normals)
                normals = vertices;
            end
            normals = normals./(sqrt(sum(normals.^2,2))*[1 1 1]);
            area1 = geometricTools.getSurfaceArea(vertices,faces);
            area2 = geometricTools.getSurfaceArea(vertices+normals,faces);
            if area2 < area1% && normalsIn
                faces = fliplr(faces);
                h = figure('visible','off');h2 = patch('vertices',vertices,'faces',fliplr(faces));
                normals = get(h2,'vertexnormals');close(h);
                normals = normals./(sqrt(sum(normals.^2,2))*[1 1 1]);
            end
            if normalsIn
                faces = fliplr(faces);
                h = figure('visible','off');h2 = patch('vertices',vertices,'faces',fliplr(faces));
                normals = get(h2,'vertexnormals');close(h);
                normals = normals./(sqrt(sum(normals.^2,2))*[1 1 1]);
            end
        end
        %%
        function [area,areas] = getSurfaceArea(vertices,faces)
            x1= vertices(faces(:,1),1);
            y1= vertices(faces(:,1),2);
            z1= vertices(faces(:,1),3);
            x2= vertices(faces(:,2),1);
            y2= vertices(faces(:,2),2);
            z2= vertices(faces(:,2),3);
            x3= vertices(faces(:,3),1);
            y3= vertices(faces(:,3),2);
            z3= vertices(faces(:,3),3);
            area = sqrt(((y2-y1).*(z3-z1)-(y3-y1).*(z2-z1)).^2+((z2-z1).*(x3-x1)-(z3-z1).*(x2-x1)).^2+...
                ((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)).^2)/2;
            areas = area;
            area = sum(area);
        end
        function L = getSurfaceLaplacian1(vertices,faces)
            % LAPLACES Calculates a Discrete Surface Laplacian Matrix
            %          for a triangulated surface
            
            Nv = size(vertices,1);
            Nf = size(faces,1);
            L = spalloc(Nv,Nv,3*Nf);
            %L = speye(Nvtx);
            for fi=1:Nf
                for k=1:3
                    kk = mod(k,3)+1;
                    L(faces(fi,k),faces(fi,kk)) = sqrt( sum((vertices(faces(fi,k),:)-vertices(faces(fi,kk),:)).^2,2));
                end
            end
            L(L>0) = 1./L(L>0);
            L = (L+L')/2;
            L = L - spdiags(sum(L,2),0,Nv,Nv);
        end
        function L = getSurfaceLaplacian(vertices,faces)
            % LAPLACES Calculates a Discrete Surface Laplacian Matrix
            %          for a triangulated surface
            %
            % Reference:
            % [1] Huiskamp, G., 1991, Difference formulas for the surface laplacian
            %     on a triangulated surface, Journal of Computational Physics 95,
            %     477-496.
            %
            % Nelson Trujillo Barreto
            % Pedro antonio Valdes Hernandez
            % Cuban Neuroscience Center
            
            Cortex.vertices = vertices;
            Cortex.faces = faces;
            vtx = Cortex.vertices;
            tri = Cortex.faces;
            [nei,nei_tri] = geometricTools.get_neis(Cortex);
            Nvtx = size(vtx,1);
            L = speye(Nvtx);
            for j=1:Nvtx,
                nei_tri_j = tri(nei_tri{j},:);
                nei_j = nei{j};
                PHI_jk = [];
                rj = vtx(j,:);
                Nj = length(nei_j);
                for k=1:Nj,
                    
                    [indi,indj]=find(nei_tri_j==nei_j(k)); %#ok
                    nei_tri_jk = nei_tri_j(indi,:);
                    if size(nei_tri_jk,1) < 2,
                        break;
                    end
                    nei_k_lr = setxor(nei_tri_jk(1,:),nei_tri_jk(2,:));
                    
                    rk2 = sum((rj-vtx(nei_j(k),:)).^2);
                    rl2 = sum((rj-vtx(nei_k_lr(1),:)).^2);
                    rr2 = sum((rj-vtx(nei_k_lr(2),:)).^2);
                    rkl2 = sum((vtx(nei_k_lr(1),:)-vtx(nei_j(k),:)).^2);
                    rkr2 = sum((vtx(nei_k_lr(2),:)-vtx(nei_j(k),:)).^2);
                    
                    cos_phi_kl = (rk2+rl2-rkl2)./sqrt(rk2.*rl2)./2;
                    cos_phi_kr = (rk2+rr2-rkr2)./sqrt(rk2.*rr2)./2;
                    sin_phi_kl = sqrt(1-cos_phi_kl.^2);
                    sin_phi_kr = sqrt(1-cos_phi_kr.^2);
                    
                    PHI_jk = [PHI_jk (1-cos_phi_kl)./(sin_phi_kl+eps)+(1-cos_phi_kr)./(sin_phi_kr+eps)]; %#ok
                end
                if ~isempty(PHI_jk)
                    rjk = sqrt(sum((vtx(nei_j,:)-repmat(rj,Nj,1)).^2,2));
                    rj_bar = mean(rjk);
                    
                    theta_jk = 4*PHI_jk'./(rj_bar*sum(PHI_jk)*rjk);
                    L(j,nei_j) = theta_jk'; %#ok
                else
                    L(j,nei_j) = 0;     %#ok
                    L(j,j) = 1;         %#ok
                end
            end;
            
            L = L - speye(Nvtx,Nvtx);
            L = L - spdiags(sum(L,2),0,Nvtx,Nvtx);
            d = diag(L);
            ind = find(d==0);
            if ~isempty(ind), for it=1:length(ind), L(ind(it),ind(it)) = 1;end;end
            %th = prctile(nonzeros(L),[0.1 99.9]);
            %L(L>th(2)) = 0;
            %L(L<th(1)) = 0;
        end
        %%
        function [nei,nei_tri] = get_neis(P)
            % helper function for getSurfaceLaplacian
            % Nelson Trujillo Barreto
            % Pedro antonio Valdes Hernandez
            % Cuban Neuroscience Center
            n = size(P.vertices,1);
            nei_tri = cell(n,1);
            nei = cell(n,1);
            for i = 1:n
                [r,c] = find(P.faces == i); %#ok
                nei_tri{i} = r;
                nei{i} = setdiff(unique(P.faces(r,:)),i);
            end
        end
        function [Nei_faces,Nei_vertices] = get_neis1(P)
            n = size(P.vertices,1);
            Nei_faces = cell(n,1);
            Nei_vertices = cell(n,1);
            % hbar = waitbar(0,'calculating neigs...');
            for i = 1:n
                [r,c] = find(P.faces == i); %#ok
                tmp = P.faces(r,:)';
                tmp(tmp == i) = [];
                m = length(tmp)/2;
                Nei_faces{i} = reshape(tmp,2,m)';
                for j = 1:m
                    Nei_vertices{i}{j} = P.vertices(Nei_faces{i}(j,:),:);
                end
                waitbar(i/n,hbar);
            end
            close(hbar);
        end
       %%
        function [nVertices,nFaces] = openSurface(vertices,faces,rmIndices)
            nVertices = vertices;
            vertices(rmIndices,:) = [];
            [~,rm_1] = ismember(faces(:,1),rmIndices);
            [~,rm_2] = ismember(faces(:,2),rmIndices);
            [~,rm_3] = ismember(faces(:,3),rmIndices);
            rm_faces = rm_1 | rm_2 | rm_3;
            faces(rm_faces,:) = [];
            [~,J] = ismember(nVertices,vertices,'rows');
            nFaces = J(faces);
            nVertices = vertices;
        end
        %%
        function [nVertices,nFaces] = getSurfaceROI(vertices,faces,roiIndices)
            rmIndices = setdiff(1:size(vertices,1),roiIndices);
            [nVertices,nFaces] = geometricTools.openSurface(vertices,faces,rmIndices);
        end
        %%
        function yi = interpOnSurface(vertices,faces,elec,y,method)
            if nargin < 5, method = 'spline';end
            switch method
                case 'ridge'
                    yi = geometricTools.ridgeInterpolation(vertices,faces,elec,y);
                case 'linear'
                    W = geometricTools.localGaussianInterpolator(elec,vertices,32);
                    yi = W*y;
                case 'spline'
                    yi = geometricTools.spSplineInterpolator(elec,y,vertices);
                otherwise
                    yi = geometricTools.spSplineInterpolator(elec,y,vertices);
            end
        end
        %%
        function atlas = labelSurface(Surf,imgAtlasfile, txtAtlasLabel,maxColorValue)
            if nargin < 4, maxColorValue = 90;end
            % Atlas
            v =spm_vol(imgAtlasfile); % atlas
            A = spm_read_vols(v);
            A(A>maxColorValue) = 0;
            indNonZero = A(:)~=0;
            colorTable = A(indNonZero);
            [x,y,z] = ndgrid(1:v.dim(1),1:v.dim(2),1:v.dim(3));
            M = v.mat;
            X = [x(:) y(:) z(:) ones(numel(x),1)]*M';
            X = X(indNonZero,1:3);           
            clear x y z
            F = TriScatteredInterp(X,colorTable,'nearest');
            n = size(Surf.vertices,1);
            labelsValue = F(Surf.vertices);
            colorTable = labelsValue;
            hwait = waitbar(0,'Atlas correction...');
            for it=1:n
                neigInd = any(Surf.faces == it,2);
                vertexInedex = Surf.faces(neigInd,:);
                vertexInedex = vertexInedex(:);
                [y,x] = hist(labelsValue(vertexInedex));
                [~,loc] = max(y);
                [~,loc] = min(abs(labelsValue(vertexInedex) - x(loc)));
                labelsValue(it) = colorTable(vertexInedex(loc));
                waitbar(it/n,hwait);
            end
            waitbar(1,hwait);
            close(hwait);
            atlas.colorTable = labelsValue;
            atlas.label = textfile2cell(txtAtlasLabel);
            atlas.label = atlas.label(1:max(atlas.colorTable));
            for it=1:length(atlas.label)
                ind = find(atlas.label{it} == ' ');
                atlas.label{it} = atlas.label{it}(ind(1)+1:ind(end)-1);
            end
        end
        %%
        function X = projectOntoUnitarySphere(X)
            [~,X(:,1),X(:,2),X(:,3)] = geometricTools.projectOnSphere(X(:,1),X(:,2),X(:,3));
            [azimuth,elevation,r] = cart2sph(X(:,1),X(:,2),X(:,3));
            [X(:,1),X(:,2),X(:,3)] = sph2cart(azimuth,elevation,elevation*0+1);
        end
        %%
        function [Yi,W] = spSplineInterpolator(X,Y,Xi,plotFlag)
            % Computes the spherical spline interpolator based on Perrin, F., 
            % Pernier, J., Bertrand, O., Echallier, J.F. (1990). Corrigenda 
            % EEG 02274. Electroencephalography and Clinical Neurophysiology, 76, 565.
            
            if nargin < 4, plotFlag = false;end
            %X0 = mean(X);
            %X  = bsxfun(@minus,X,X0);
            %Xi  = bsxfun(@minus,Xi,X0);
            X  = geometricTools.projectOntoUnitarySphere(X);
            Xi = geometricTools.projectOntoUnitarySphere(Xi);
            %X  = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
            %Xi = bsxfun(@rdivide,Xi,sqrt(sum(Xi.^2,2)));
            
            %-- 
            M = size(X,1);
            One = ones(size(Xi,1),1);
            %--
            
            % Solving eq. 4 of Perrin et al. (1989)
            COS_X  = geometricTools.cosines(X,X);
            COS_Xi = geometricTools.cosines(Xi,X);
            
            % Solving eq. 3 of Perrin et al. (1989)
            Gx  = geometricTools.sphericalSpline(COS_X);
            Gxi = geometricTools.sphericalSpline(COS_Xi);
            
            % Solving eq. 2 Perrin et al. (1989)
            [C,~,~,T] = ridgeGCV([Y;0],[Gx ones(M,1);ones(1,M) 0],eye(M+1));
                        
            % Interpolating with the spherical harmonics
            Yi = [Gxi One]* C;
            
            W = [Gxi One]* T(:,1:end-1);
            
            % Plot the input & projected electrode positions on a sphere
            if plotFlag
                geometricTools.plot_on_sphere(X,Y,Xi,Yi);
            end
        end
        %%
        function Gx = sphericalSpline(x)
            % sphericalSpline solves eq. 3 of Perrin et al. (1989)
            % g(COS) = 1/4pi * sum[n=1:inf] (( (2*n+1)/( n^m * (n+1)^m ) ) * Pn(COS));
            
            m = 4;
            N = 16;    % gives accuracy of 10^-6
            
            P = cat(3, ones(size(x)), x);
            Gx = 3 / 2 ^ m * P(:, :, 2);
            for n = 2:N
                P(:, :, 3) = ((2 * n - 1) * x .* P(:, :, 2) - (n - 1) * P(:, :, 1)) / n;
                P = P(:,:,[2 3 1]);
                Gx = Gx + (2 * n + 1) / (n ^ m * (n + 1) ^ m) * P(:, :, 2);
            end
            Gx = Gx / (4 * pi);
        end
        %%
        function [r,x,y,z] = projectOnSphere(X,Y,Z,xo,yo,zo)
            % projectOnSphere - calculates projections of xyz positions
            % onto the unitary sphere
            %
            % Usage: [r,x,y,z] = projectOnSphere(X,Y,Z,xo,yo,zo,plotFlag)
            %
            % Notes:    The general formula for a sphere, with radius r is given by:
            %
            %           (x - xo)^2  +  (y - yo)^2  +  (z - zo)^2  =  r^2
            %
            %           This function takes arguments for cartesian co-ordinates
            %           of X,Y,Z (assume Z > 0) and the center of the sphere (xo,yo,zo).
            %           If (xo,yo,zo) is not provided a cnter at (0,0,0) is assumed.
            %
            %           Returned values are the fitted radius 'r' (constant)
            %           and the (x,y,z) Cartesian coordinates of the projected points
            %
            %
            % $Revision: 1.3 $ $Date: 2005/07/12 22:16:48 $
            % Licence:  GNU GPL, no express or implied warranties
            % History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
            %                    adapted from elec_fit_sph
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % initialise centroid, unless input parameters defined
            if nargin < 4, xo = 0;end
            if nargin < 5, yo = 0;end
            if nargin < 6, zo = 0;end
            
            % Initialise r0 as a rough guess at the sphere radius
            rX = (max(X) - min(X)) / 2;
            rY = (max(Y) - min(Y)) / 2;
            rZ =  max(Z) - zo;
            r0 = mean([ rX rY rZ ]);
            
            % perform least squares estimate of spherical radius (r)
            options = optimset('fminsearch');
            r = fminsearch(@geometricTools.fit2sphere,r0, options, X, Y, Z, xo, yo, zo);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the projection point of X,Y,Z to the fitted sphere radius r
            
            % Convert Cartesian X,Y,Z to spherical (radians)
            theta = atan2( (Y-yo), (X-xo) );
            phi = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
            % do not recalc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
            
            %   Recalculate X,Y,Z for constant r, given theta & phi.
            R = ones(size(phi)) * r;
            x = R .* sin(phi) .* cos(theta);
            y = R .* sin(phi) .* sin(theta);
            z = R .* cos(phi);
        end 
    end
    methods(Static,Hidden=true)
        %%
        function f = fit2sphere(r, X, Y, Z, xo, yo, zo)
            S = (X-xo).^2  +  (Y-yo).^2  +  (Z-zo).^2  -  r^2;
            f = sum( S.^2 );
        end
        %%
        function Cos = cosines(A,B)
            Na = size(A,1);
            Nb = size(B,1);
            One1 = ones(1,Nb);
            One2 = ones(Na,1);
            Xe = A(:,1)*One1;
            Ye = A(:,2)*One1;
            Ze = A(:,3)*One1;
            Xf = One2*B(:,1)';
            Yf = One2*B(:,2)';
            Zf = One2*B(:,3)';
            Cos = (Xe-Xf).^2 + (Ye-Yf).^2 + (Ze-Zf).^2;
            Cos = 1-Cos/2;
            Cos(Cos > 1) = 1-eps;
            Cos(Cos < -1) = -1+eps;
        end
        function plot_on_sphere(X,Y,Xi,Yi)
            [~,X(:,1), X(:,2), X(:,3)]  = geometricTools.projectOnSphere(X(:,1), X(:,2), X(:,3));
            [~,Xi(:,1),Xi(:,2),Xi(:,3)] = geometricTools.projectOnSphere(Xi(:,1),Xi(:,2),Xi(:,3));
            X  = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
            Xi = bsxfun(@rdivide,Xi,sqrt(sum(Xi.^2,2)));
            
            Xt = [X;Xi];
            Yt = [Y;Yi];
            Xt = bsxfun(@rdivide,Xt,sqrt(sum(Xt.^2,2)));
            Xi = bsxfun(@rdivide,Xi,sqrt(sum(Xi.^2,2)));
            Ne = size(X,1);
            Nf = 72;
            [Xs,Ys,Zs]=sphere(Nf);
            Xsp = [Xs(:) Ys(:) Zs(:)];
            Fsp = geometricTools.localGaussianInterpolator(Xt,Xsp,0.2);
            
            
            %[J,lambdaOpt,~,iFsp] = ridgeGCV(Yt,Fsp',eye(size(Xsp,1)),100,1);
            %J = iFsp*Yt;
            Ysp = Fsp*Yt;
            figure('NumberTitle','off','Name','Electrode Placements');
            set(gca,'Projection','perspective','DataAspectRatio',[1 1 1]); hold on
            %plot3(x,y,z,'b.');
            plot3(X(:,1),X(:,2),X(:,3),'ro');
            plot3(Xi(:,1),Xi(:,2),Xi(:,3),'k.')
            legend('input xyz','projected head','Location','BestOutside');
            surf(Xs,Ys,Zs,reshape(Ysp,[Nf Nf]+1),'specularstrength',0.1,'facealpha',0.9,'linestyle','none');
            camlight
            camlight headlight
            view(2); rotate3d;
            axis vis3d
        end
    end
end