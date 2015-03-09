load ../cognionicsHeadModel
load ../cognionicssLORETA
load(['../' cognionicsHeadModel.surfaces])
K_full = zeros(size(cognionicsHeadModel.channelSpace,1),size(surfData(3).vertices,1)*3);
K_full(:,[ind;ind+n;ind+2*n]) = K;


pf.vertices = surfData(3).vertices;
pf.g_matrix = K_full;
pf.sigma = 1; % not sure
pf.cov_matrix = eye(size(pf.g_matrix,1));
pf.ssp = eye(size(pf.g_matrix,1));
pf.matrix_pw = eye(size(pf.g_matrix,1)); % ????
pf.evol = eye(size(pf.g_matrix,2));
pf.np = 5000;
pf.sigma_par = 1; % not sure
pf.zero_time = 0; % not sure
pf.final_time = 100; % not sure
pf.freq = 10; % not sure
pf.time_interval = pf.zero_time/1000:1/pf.freq:pf.final_time/1000;
pf.t1 = 1;
pf.t2 = 2;
pf.data =    1.0e-05 *[0.7396
                    0.8303
                    0.8626
                    0.8964
                    0.8575
                    0.7007
                    0.7829
                    0.7032
                    0.8599
                    0.8744
                    0.8801
                    0.8385
                    0.7953
                    0.6382
                    0.3786
                    0.5738
                    0.6679
                    0.6924
                    0.6889
                    0.6762
                    0.6077
                    0.5301
                    0.3381
                    0.0937
                    0.2958
                    0.3368
                    0.4347
                    0.3707
                    0.3621
                    0.2994
                    0.1839
                    0.0228
                    -0.2713
                    -0.0970
                    -0.0713
                    -0.0494
                    -0.0726
                    -0.0433
                    -0.0637
                    -0.1049
                    -0.0686
                    -0.4382
                    -0.4907
                    -0.5011
                    -0.5861
                    -0.5223
                    -0.5059
                    -0.4683
                    -0.4012
                    -0.3553
                    -0.7029
                    -0.7976
                    -0.9290
                    -0.9233
                    -0.9263
                    -0.7331
                    -0.6287
                    -0.7181
                    -0.8406
                    -0.9074
                    -0.8609
                    -0.8866
                    -0.8213
                    -0.6296]';
pf.data(2,:) = 0.9*pf.data(1,:);