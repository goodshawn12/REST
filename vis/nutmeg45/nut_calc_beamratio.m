function [beam,selnum] = nut_calc_beamratio(beam)
% Calculates ratio of values in beam.s
% Calls a GUI to let user select which ratio.
%  beam = nut_calc_beamratio(beam) 

global stats

if ~isfield(stats,'ratio')
    gh=nut_select_beamratio;
    uiwait(gh)
end

if stats.noise
    noise=beam.s{3};
else
    noise=0;
end

switch stats.ratio
    case 1  % 'F-ratio (dB)'
        s{1} = 10*log10(abs((beam.s{1}-noise)./(beam.s{2}-noise)));
    case 2  % 'F-ratio (raw)'
        s{1} = (beam/s{1}-noise)./(beam.s{2}-noise);
    case 3  % 'CTF Pseudo-F'
        F = (beam.s{1}-noise)./(beam.s{2}-noise);
        selectpos = (F>=1);
        selectneg = (F<1 & F>0);
        % selectzero = find(rivets.s==0);
        s{1}=zeros(size(F));
        s{1}(selectpos) = F(selectpos) - 1;
        s{1}(selectneg) = 1-(1./F(selectneg));
    case 4  % 't-test (difference)'
        s{1} = (beam.s{1}-noise) - (beam.s{2}-noise);
    case 5  % change
        s{1} = 100*(beam.s{1} - beam.s{2})./(beam.s{2}-noise);
    case 6  % Wilcoxon Z score
        % s{1} = beam.z;
        %% note that spatnorm files currently place beam.z into beam.s{1} for convenience
        s{1} = beam.s{1};
    case 7   % Active
        s{1} = beam.s{1}-noise;
    case 8   % Control
        s{1} = beam.s{2}-noise;
    case 9   % Imaginary Coherence
        s{1} = beam.s{1};
    case 10  % Magnitude Squared Coherence
        s{1} = beam.s{2};
    case 11
        s{1} = beam.s{3};
end
beam.s=s;

labels = {'F-ratio (dB)' 'F-ratio (raw)' 'CTF Pseudo-F' '\Delta Power' 'Power Change (%)' ...
    'Wilcoxon Z score' 'Active' 'Control' 'Imaginary Coherence' 'Magnitude Squared Coherence' 'Real Coherence'};
beam.sinfo = labels(stats.ratio);

if nargout>1, selnum=stats.ratio; end
