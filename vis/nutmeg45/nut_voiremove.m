function voi=nut_voiremove(voi,flag,number)
% voi=nut_voiremove(voi,flag,number)
%
% Removes components from voi structure
% flag      'c'=condition, 's'=subject, 'voi'=VOI
% number    flag number(s) to be removed.
%

switch flag
    case {'c' 'cond' 'condnr'}
        f=find(ismember(voi.condnr,number));
        voi.subjnr(f)=[];
        voi.condnr(f)=[];
        voi.filenames(f)=[];
        voi.pathnames(f)=[];
        for cc=1:length(voi.coords)
            voi.coords(cc).MRI(f,:)=[];
            voi.coords(cc).MEG(f,:)=[];
            voi.coords(cc).index(:,f)=[];
        end
   case {'s' 'subj' 'subjnr'}
        f=ismember(unique(voi.subjnr),number);
        voi.subjectmaxtf(f)=[];
        f=find(ismember(voi.subjnr,number));
        voi.subjnr(f)=[];
        voi.condnr(f)=[];
        voi.filenames(f)=[];
        voi.pathnames(f)=[];
        for cc=1:length(voi.coords)
            voi.coords(cc).MRI(f,:)=[];
            voi.coords(cc).MEG(f,:)=[];
            voi.coords(cc).index(:,f)=[];
        end
    case {'coord' 'voi'}
        voi.coords(number)=[];
end