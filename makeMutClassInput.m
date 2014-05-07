addpath('/home/kchang3/data/ill-snv/');
originalGenomes=importdata('prim-originalGenomes.txt')
subtypes=importdata('subtypes.txt')
types=importdata('types.txt')
sampleNames=importdata('prim-sampleNames.txt')
cancerType='ova triplets' 
%strandTypes=load('strandTypes.txt')

save('/home/kchang3/data/ill-snv/prim.mat','originalGenomes','subtypes','types','sampleNames','cancerType')
originalGenomes=importdata('recur-originalGenomes.txt')
sampleNames=importdata('recur-sampleNames.txt')
save('/home/kchang3/data/ill-snv/recur.mat','originalGenomes','subtypes','types','sampleNames','cancerType')
originalGenomes=importdata('originalGenomes.txt')
sampleNames=importdata('sampleNames.txt')
save('/home/kchang3/data/ill-snv/all.mat','originalGenomes','subtypes','types','sampleNames','cancerType')

quit
