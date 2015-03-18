function Fsummary=nut_find_same_voxels(VOIvoxels, beam, timewin, freqbin, k, d)
%function voxel_mask=nut_find_same_voxels(VOIvoxels, beam, timewin, freqbin)
% VOIvoxels = variable stored in voxel file made using Nutmeg VOI tool
% beam = variable stored in s_beam... file
% timewin = number of time window used as in time-freq plot (1 will be
%       earliest time window, 2 next etc.)
% freqbin = number of frequency bin used as in time-freq plot (1 will be
%       lowest frequency bin, 2 next etc.)
% k is the same size as beam.voxel. The index # of k corresponds to the
% voxel number in beam.voxel. The value of k(i) corresponds to a voxel in
% the VOI with distance d(i) from the kth voxel in the beam.voxel. Does
% that make sense?
%
% written by Anne Findlay, 2008

voxelsize = beam.voxelsize;

%%Uncomment this section if you are not using k & d as inputs
%for ii = 1:3
%    voxels(:,ii) = voxelsize(ii)*round(VOIvoxels(:,ii)/voxelsize(ii));
%end

%uniquevoxels = unique(voxels,'rows');

%[k d] = dsearchn(uniquevoxels,beam.voxels);

for timewin = 1:size(beam.s{1},2) %iterate over time windows, earliest first
    for freqbin = 1:size(beam.s{1},3) %iterate over frequency bins, lowest first

        for ii = 1:size(k,1)
            if d(ii) <= 10
                voxel_number(ii) = k(ii); %k(ii) is the index number of the VOI
            elseif d(ii) > 10
                voxel_number(ii) = 0;
            end
        end

        voxel_mask(1:size(k,1)) = 0;

        for ii = 1:size(k,1)
            if voxel_number(ii) ~= 0
                voxel_mask(ii) = 1;
            end
        end

        voxel_mask = voxel_mask';

        masked_beam = voxel_mask.*beam.s{1}(:,timewin,freqbin);

        %To test your output....
        %beam.s{1}(:,1,1)=voxel_mask.*beam.s{1}(:,1,1);
        %save testbeam.mat beam

        Fsummary(timewin,freqbin).averageF = sum(masked_beam)/sum(voxel_mask)
        Fsummary(timewin,freqbin).maxF=max(masked_beam)
        Fsummary(timewin,freqbin).minF=min(masked_beam)

    end
end