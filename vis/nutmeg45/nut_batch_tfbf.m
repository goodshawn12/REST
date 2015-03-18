% nut_savefilternuts('allinfo.mat',8,15,'firlsbp200cn.mat')
% nut_constructTFBFparams(0,6700,[100 150 200 300],50,7400,7500);
% nut_constructTFBFparams([50 100 150 250],[19600 19550 19500 19400],[100 150 200 300],50,-300,-200);
% nut_calltfcov('allinfo.mat','300ms',8,15,'firlsbp200cn.mat')
% nut_liposession('allinfo.mat');
% nut_calltfrun('allinfo.mat','600ms',15,30,'firlsbp200cn.mat')

nd=3;
dir{1}='c1_450_16700';
dir{2}='c2_450_16700';
dir{3}='c3_450_16700';
% freq=[8 15; 15 30; 30 45; 55 95; 105 145; 155 195];
% freq=[15 30; 30 45; 55 95; 105 145; 155 195];
% freq=[155 195];
freq=[0 2;2 4;4 8];
etime=16700;
% windur=[450 300 200 150 150 150]; %must be same length as size(freq,1)
% windur=[300 200 150 150 150]; %must be same length as size(freq,1)
% windur=[150]; %must be same length as size(freq,1)
windur=[1800 1200 600];
for jj=1:size(freq,1);
    for ii=1:nd
        aa=load(['../' dir{ii} '/allinfo_c' num2str(ii) 'unave_firlsbp200cn_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_con' num2str(etime-windur(jj)) 'to' num2str(etime) 'ms.mat']);
        if ii==1
            R=zeros(size(aa.R));
            Rcon=zeros(size(aa.Rcon));
        end
        R=R+aa.R;
        Rcon=Rcon+aa.Rcon;
        params=aa.params;
        clear aa
    end
    save(['allinfo_c123unave_firlsbp200cn_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_con' num2str(etime-windur(jj)) 'to' num2str(etime) 'ms.mat'],'R','Rcon','params');
    clear R Rcon
    nut_calltfrun('allinfo_c123unave.mat',[num2str(windur(jj)) 'ms'],freq(jj,1),freq(jj,2),'firlsbp200cn.mat');
end

windur=[450];freq=[52 98];jj=1;covfile_usechar{1}=5:30;
nut_calltfrun('all_c123.mat',[num2str(windur(jj)) 'ms'],freq(jj,1),freq(jj,2),'firlsbp200cn.mat',covfile_usechar{jj});



% tfcov('allinfo.mat','300ms.mat',8,15,'firlsbp200cn.mat',{'SAM'});
% tfcov('allinfo.mat','200ms.mat',15,30,'firlsbp200cn.mat',{'SAM'});
% tfcov('allinfo.mat','150ms.mat',30,45,'firlsbp200cn.mat',{'SAM'});
% tfcov('allinfo.mat','100ms.mat',55,95,'firlsbp200cn.mat',{'SAM'});
% tfcov('allinfo.mat','100ms.mat',105,145,'firlsbp200cn.mat',{'SAM'});
% tfcov('allinfo.mat','100ms.mat',155,195,'firlsbp200cn.mat',{'SAM'});
% 
% tfcov('allinfo.mat','300ms.mat',8,15,'firlsbp200cn.mat');
% tfcov('allinfo.mat','200ms.mat',15,30,'firlsbp200cn.mat');
% tfcov('allinfo.mat','150ms.mat',30,45,'firlsbp200cn.mat');
% tfcov('allinfo.mat','100ms.mat',55,95,'firlsbp200cn.mat');
% tfcov('allinfo.mat','100ms.mat',105,145,'firlsbp200cn.mat');
% tfcov('allinfo.mat','100ms.mat',155,195,'firlsbp200cn.mat');
% 
% tfcov('allinfo_c3unave.mat','300ms.mat',8,15,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave.mat','200ms.mat',15,30,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave.mat','150ms.mat',30,45,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave.mat','100ms.mat',55,95,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave.mat','100ms.mat',105,145,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave.mat','100ms.mat',155,195,'firlsbp200cn.mat','SAM');
% 
% nut_tfbf2timef allinfo_c3unave_firlsbp200cn SAM
% 
% nut_liposession('allinfo_c3unave3gr.mat');
% 
% tfcov('allinfo_c3unave3gr.mat','300ms.mat',8,15,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','200ms.mat',15,30,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','150ms.mat',30,45,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','100ms.mat',55,95,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','100ms.mat',105,145,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','100ms.mat',155,195,'firlsbp200cn.mat','SAM');
% 
% nut_tfbf2timef allinfo_c3unave3gr_firlsbp200cn SAM
% 
% 
% tfcov('allinfo_c3unave3gr.mat','600ms.mat',8,15,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','400ms.mat',15,30,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','300ms.mat',30,45,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','200ms.mat',55,95,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','200ms.mat',105,145,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unave3gr.mat','200ms.mat',155,195,'firlsbp200cn.mat','SAM');
% 
% nut_tfbf2timef allinfo_c3unave3gr_firlsbp200cn SAM
% 
% tfcov('allinfo_c3unaveChrm1-3gr.mat','600ms.mat',8,15,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unaveChrm1-3gr.mat','400ms.mat',15,30,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unaveChrm1-3gr.mat','300ms.mat',30,45,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unaveChrm1-3gr.mat','200ms.mat',55,95,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unaveChrm1-3gr.mat','200ms.mat',105,145,'firlsbp200cn.mat','SAM');
% tfcov('allinfo_c3unaveChrm1-3gr.mat','200ms.mat',155,195,'firlsbp200cn.mat','SAM');
% 
% nut_tfbf2timef allinfo_c3unaveChrm1-3gr_firlsbp200cn SAM
% 


