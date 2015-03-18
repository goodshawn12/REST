%%
% Finds the Berg parameters using a Nelder-Mead simplex search. Algorithm
% based on Zhang's paper.
% @param r - sphere radii from outside-in [mm].
% @param cond - conductivities from outside-in [S/m]
% @see Zhang, Phys. Med. Biol., 1995
% @author Daniel D.E. Wong
%%
function [mu,lambda] = nut_findbergparams(r,cond)
maxIt = 25;

% Computations in Zhang's paper assume r and ci indexed inside-outside
rkn = flipud(fliplr(r))/1000;
ci = flipud(fliplr(cond));

M = length(r);

N = 40;
fn = zeros(N,1);
for n = 1:N
    m = 1/(2*n+1)^(M-1);
    for k = 1:M-1
        m = m * [n+(n+1)*ci(k)/ci(k+1), (n+1)*(ci(k)/ci(k+1)-1)*(rkn(end)/rkn(k))^(2*n+1); ...
            n*(ci(k)/ci(k+1)-1)*(rkn(k)/rkn(end)).^(2*n+1), (n+1)+n*ci(k)/ci(k+1)];
    end
    fn(n) = n/(n*m(2,2) + (1+n)*m(2,1));
end


options.tolX = 1e-3;
options.maxIt = 1e3;
cflag = 0;
it = 0;
warning off
while ~cflag
    it = it + 1; if it > maxIt; break; end;
    mulambda = 2*(rand([M,2,M*2+1])-0.5);
    [x,delta,cflag] = nut_simplexSearch(@optfcn,mulambda,options,rkn,fn);
end
warning on
if it > maxIt; warning('Could not find berg parameters'); end;

mu = x(:,1);
lambda = x(:,2);
lambda(1) = fn(1) - sum(lambda(2:end));


%%%%%%%%%%%%% Test Code %%%%%%%%%%%%%%%
% vol.r = r;
% vol.cond = cond;
% deg = [0:5:359]';
% sensorcoord = [r(1)*cosd(deg) r(1)*sind(deg) zeros(length(deg),1)]/1000;
% rq = [0 0.8*r(end) 0]/1000;
% lx1 = zeros(size(sensorcoord,1),1);
% for k = 1:size(x,1)
%     for s = 1:length(lx1)
%         [lx1_,ly1_,lz1_] = nut_1lyrsphleadp(rq*1000*mu(k),sensorcoord(s,:)*1000,[0 0 0],vol);
%         lx1(s) = lx1(s) + lambda(k)*lx1_;
%     end
% end
% figure, plot(lx1)
% 
% vol.r = r;
% vol.cond = cond;
% lx1_0 = zeros(size(sensorcoord,1),1);
% for s = 1:length(lx1_0)
%     [lx1_,ly1_,lz1_] = nut_sphleadp(rq*1000,sensorcoord(s,:)*1000,[0 0 0],vol);
%     lx1_0(s) = lx1_;
% end
% figure, plot(lx1_0)



function delta = optfcn(mulambda,rkn,fn)
delta = 0;
for n = 2:length(fn)
    ml = mulambda(:,2).*(mulambda(:,1).^(n-1)-mulambda(1,1)^(n-1));
    ml = sum(ml,1);
    delta = delta+((rkn(1)/rkn(min(n,length(rkn))))^(n-1)*(fn(n) - fn(1)*mulambda(1,1)^(n-1) - ml))^2;
end
