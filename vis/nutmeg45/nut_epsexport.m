function nut_epsexport(filename,figh)
% exports figures as EPS, with adobecset workaround for missing
% dashes/minus signs

if(~exist('figh','var'))
    figh=gcf;
end

print(figh,'-depsc2','-adobecset',filename)
