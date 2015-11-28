clinical_file = 'clinical2_brca.csv';

fid = fopen(clinical_file);
%                    1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 
data = textscan(fid,'%s %d %f %f %f %f %s %s %s %s %s %s %s %s', 'delimiter', ',');
ppheno.sample_id = data{1};
ppheno.survival = data{2};
ppheno.days_to_death = data{3};
ppheno.days_to_last = data{4};
ppheno.days_survival = data{5};
ppheno.diag_age = data{6};
ppheno.subdivision = data{7};
ppheno.hist_type = data{8};
ppheno.pathologic_stage = data{9};
ppheno.cancer_status = data{10};
ppheno.estrogen_status = data{11};
ppheno.progesterone_status = data{12};
ppheno.margin_status = data{13};
ppheno.menopause_status = data{14};
fclose(fid);

clear data;
clear fid;
clear filtered_id;
clear basepath;
clear clinical_file;
clear ans;

