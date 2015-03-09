function[]=hades_write_BS
global struct_ris_4BS;
global pf;
global ris;
global gui_bs;

% matrix n_vertici X n_time_frame
ris_for_BS = zeros(size(pf.vertices,1), length(ris.time_interval));
% allocate in the first time point all reconstructed sources
for i=1:size(ris.dip,1)
    i_tf = ris.dip(i,end-1);
    i_vertix = ris.dip(i,end);
    ris_for_BS(i_vertix, 1) = ris_for_BS(i_vertix, 1) + 1;
    ris_for_BS(i_vertix, i_tf) = 1;
end
ris_for_BS(:,1) = ris_for_BS(:,1)/max(ris_for_BS(:,1));

% create struct to import ris data in BS
% these variables have to be written by HADES
struct_ris_4BS.ChannelFlag = ones(size(ris.data,2),1);
struct_ris_4BS.Comment = 'PF results';
struct_ris_4BS.Function = 'pf';
struct_ris_4BS.GoodChannel = 1:size(ris.data,2);
struct_ris_4BS.ImageGridAmp = ris_for_BS;
struct_ris_4BS.ImageGridTime = ris.time_interval;
struct_ris_4BS.ImagingKernel = [];
struct_ris_4BS.Time = 1:length(struct_ris_4BS.ImageGridTime);
struct_ris_4BS.Whitener = pf.matrix_pw;

% these variables have to be set in BRAINSTORM
struct_ris_4BS.DataFile = get(gui_bs.str_DataFile,'string');
struct_ris_4BS.HeadModelFile = get(gui_bs.str_HeadModelFile,'string');
struct_ris_4BS.SurfaceFile = get(gui_bs.str_SurfaceFile,'string');
if strcmp(struct_ris_4BS.DataFile,'')==0 && strcmp(struct_ris_4BS.HeadModelFile,'')==0 && strcmp(struct_ris_4BS.SurfaceFile,'')==0
    if exist(ris.results_dir)==0
        start_path=pwd;
        directory_name = uigetdir(start_path,'Select the directory where save the BrainStorm files');
    else
        directory_name = ris.results_dir;
    end
    filename_bs = fullfile(directory_name,strcat(ris.subject,'_bs.mat'));
    save(filename_bs,'struct_ris_4BS','-mat');
    close(findobj('Name','Export data to BrainStorm','-and','-depth',6));
else
    msgbox('Some fields have not been filled','Error','error');
end
    

return