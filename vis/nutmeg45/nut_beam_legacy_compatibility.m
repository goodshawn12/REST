function beam=nut_beam_legacy_compatibility(beam)
% function beam=nut_beam_legacy_compatibility(beam)
%
% copied from nut_results_viewer.m on 25 Aug, 2010

if isfield(beam,'beam')
    beam=beam.beam;
end
if ~isfield(beam,'s')
    if(isfield(beam,'tf'))  
        beam.s{1}=beam.tf;
        beam = rmfield(beam,'tf');
    elseif isfield(beam,'s_th') % compatible with nut_activation_viewer
        beam.timepts = beam.timewindow;
        beam = rmfield(beam,'timewindow');
        beam.s{1}(:,:,1,1)=beam.s_th;
        if ~all(beam.s_ph(:)==0)
            beam.s{1}(:,:,1,2)=beam.s_ph;
        end
        beam = rmfield(beam,'s_th');
        beam = rmfield(beam,'s_ph');
        if isfield(beam,'s_z')
            if ~all(beam.s_z(:)==0)
                beam.s{1}(:,:,1,3)=beam.s_z;
            end
            beam = rmfield(beam,'s_z');
        end    
    elseif isfield(beam,'sa')
        beam.s(:,1) =  beam.sa;
        beam=rmfield(beam,'sa');
        if isfield(beam,'sc')
            beam.s(:,2) =  beam.sc;
            beam=rmfield(beam,'sc');
        end
        if isfield(beam,'n')
            beam.s(:,3) =  beam.n;
            beam=rmfield(beam,'n');
        end
    end

    if isfield(beam,'saeach')
        for ii=1:size(beam.saeach,2)
        beam.s{:,ii*2+3}=beam.saeach{ii};
        beam.s{:,ii+2+4}=beam.sceach{ii};
        end
        beam=rmfield(beam,'saeach');
        beam=rmfield(beam,'sceach');
    end
end
if isfield(beam,'ttestFDR')
    beam.ttest=beam.ttestFDR;
    beam=rmfield(beam,'ttestFDR');
end

    