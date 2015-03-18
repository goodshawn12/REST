function position=nut_vol2cloud6(thld,vol)
%function POSITION=NUT_VOL2CLOUD6(THLD,VOL)

h = waitbar(0,'Please wait, finding scalp surface...');

vol=vol(:,:,:)>thld;
%do the segmentation
[L,NUM] = bwlabeln(vol);
objects=hist(reshape(L,size(L,1)*size(L,2)*size(L,3),1),NUM+1);
%get rid of spackles and wrap around
label_of_interest=find(objects>10000); 
if objects(1)<objects(2)
	label=label_of_interest(1);
else
	label=label_of_interest(2);
end
L=(L(:,:,:)==label-1) ;     
position=zeros(0,3);
%get only the first and the last non zero value from each row of each slice
for k=1:size(vol,3) 
	[crap,p]=max(L(:,:,k));
	no_of_rows=size(vol,1);
	p=p(:);
	p_index=find(p~=1);
	[crap,p1]=max(flipud(L(:,:,k)));
	p1_real=no_of_rows+1-p1;
	p1_real=p1_real(:);
	p1_real_index=find(p1_real~=no_of_rows);
	position1=[[p(p_index),p_index,(ones(length(p_index),1)*k)];[p1_real(p1_real_index),p1_real_index,(ones(length(p1_real_index),1))*k]];
	position = [position;position1]; 
	waitbar(k/(size(vol,3)),h);
end

close(h);
