function [absdataout,angdataout]=nut_abshilbert(datain)
% [absdataout,angdataout]=nut_abshilbert(datain)
%
% datain size: time x channels x trials

absdataout = zeros(size(datain));
angdataout = zeros(size(datain));
hildata = zeros(size(datain,1),size(datain,2));

for kk=1:size(datain,3)
    hildata=hilbert(datain(:,:,kk));
    absdataout(:,:,kk)=abs(hildata);
    if nargout>1
        angdataout(:,:,kk)=angle(hildata);
    end
end

