function voi = nut_voinew(subjnr,condnr,plist,flist)
% NUT_VOINEW creates new pointerfile structure.
%
%  voi = nut_voinew(subjnr,condnr,plist,flist)
%
% SUBJNR    number of each subject (e.g.: [1:4 7:10])
% CONDNR    number of each condition
% PLIST     cell with path to data directory of each subject. For each subjects, the 
%           data of all conditions must be in the same directory.
% FLIST     cell with file name of s_beam*.mat file for each condition. All
%           subjects of a given condition must have the same file name.

if ~iscell(flist), flist={flist}; end       % if only one condition, input may come without the cell.
if length(plist)~=length(subjnr)
    error('PLIST must have same number of cells as number of subjects.')
end
if length(flist)~=length(condnr)
    error('FLIST must have same number of cells as number of conditions.')
end  

numsubj=length(subjnr);
numcond=length(condnr);

voi=struct('subjnr',[],'condnr',[],'filenames',[],'pathnames',[],'coords',[],'subjectmaxtf',[]);
for s=1:numsubj
    for c=1:numcond
        voi.subjnr(end+1)=subjnr(s);
        voi.condnr(end+1)=condnr(c);
        voi.pathnames{end+1}=[plist{s} filesep flist{c}];
        voi.filenames{end+1}=flist{c};
    end
end
voi.subjectmaxtf=zeros(1,numsubj);