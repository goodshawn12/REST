function liposession(sessionfile)

nuts=load(sessionfile);
if(isfield(nuts,'VOIvoxels'))
    nuts=rmfield(nuts,'VOIvoxels');
end
if(isfield(nuts,'meg'))
    if(isfield(nuts.meg,'data'))
        nuts.meg=rmfield(nuts.meg,'data');
    end
end

[sessionpath,sessionname,ext]=fileparts(sessionfile);
outname=[sessionname 'Lp'];
save(fullfile(sessionpath,[outname ext]),'nuts');