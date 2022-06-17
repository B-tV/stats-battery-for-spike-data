proc import datafile="C:\Mehmet_Kocak\BERD Clinic\Ben\Ready for SAS.xlsx"
dbms=xlsx out=dopa_raw replace; run;
proc import datafile="C:\Mehmet_Kocak\BERD Clinic\Ben\mehmet dataset 10-27-16.xlsx"
dbms=xlsx out=dopa_new replace; run;
proc format;
value depth 0='Not High' 1='High';
value genot 0='WT' 1='PitX'; run;
data dopa_new; set dopa_new;
depth=2-catDV2;
rename catGT=genotype; run;
proc mixed data=dopa_new plots=all covtest;
class animal genotype depth;
model logpostprerat=genotype|depth/solution ddfm=kr outpm=predictions;
repeated/type=cs subject=animal r rcorr;
run;
proc tabulate data=predictions; class genotype depth; var pred StdErrPred; format genotype genot. depth depth.;
table genotype*depth, n (pred StdErrPred)*mean*f=10.3; run;
proc tabulate data=predictions; class genotype depth; var pred StdErrPred; format genotype genot. depth depth.;
table depth*genotype, n (pred StdErrPred)*mean*f=10.3; run;
*** Covariance Structure is not very strong so we are testing the Variance Component Version ***;
proc mixed data=dopa_new plots=all covtest;
class animal genotype depth;
model logpostprerat=genotype|depth/solution ddfm=kr outpm=predictions;
repeated/type=vc subject=animal r rcorr;
run;
proc tabulate data=predictions; class Exclude genotype depth; 
format genotype genot. depth depth.;
table depth*genotype, Exclude*n; run;
proc mixed data=dopa_new plots=all; where Exclude="'Y'" and animal^="'mstrSNrfront'";
class animal genotype depth;
model logpostprerat=genotype/solution ddfm=kr outpm=predictions;
run;


proc gplot data=predictions;
symbol i=none v=circle;
plot pred*logpostprerat; run; quit;

*** Earlier Discussions ***;
data dopa; set dopa_raw;
if b='Y' then genotype=1;
else genotype=0;
depth=catDV2-1;
*if logpostprerat<-2.0 then logpostprerat=-2;
if logpostprerat>0 then response=1;
else response=0;
run;
data dopa_est; 
genotype=0; depth=0; base_ave1=0.01; output;
genotype=0; depth=1; base_ave1=0.01; output;
genotype=1; depth=0; base_ave1=0.01; output;
genotype=1; depth=1; base_ave1=0.01; output;
genotype=0; depth=0; base_ave1=1; output;
genotype=0; depth=1; base_ave1=1; output;
genotype=1; depth=0; base_ave1=1; output;
genotype=1; depth=1; base_ave1=1; output;
genotype=0; depth=0; base_ave1=0.1; output;
genotype=0; depth=1; base_ave1=0.1; output;
genotype=1; depth=0; base_ave1=0.1; output;
genotype=1; depth=1; base_ave1=0.1; output;
run; 
data dopa_est; set dopa_est; specialcases=1; run;
data dopa_all; set dopa dopa_est; run;
proc univariate data=dopa;
var logpostprerat; histogram logpostprerat; run;

proc mixed data=dopa plots=all;
class animal genotype depth;
model logpostprerat=genotype|depth|base_ave1/solution ddfm=kr;
repeated/type=vc subject=animal r rcorr; run;

proc mixed data=dopa plots=all method=reml;
class animal genotype (ref='0') depth (ref='0');
model logpostprerat=genotype|depth depth|base_ave1/solution ddfm=kr;
repeated/type=vc subject=animal r rcorr;
contrast 'Genotype' genotype 1 -1;
contrast 'Depth' depth 1 -1;
run;

proc mixed data=dopa_all;
class animal genotype (ref='0') depth (ref='0');
model logpostprerat=genotype|depth base_ave1/solution ddfm=kr outpm=estimates outp=estimates1;
repeated/type=vc subject=animal r rcorr;
run;

proc sql; select distinct genotype, depth, mean(pred) as mean_pred from estimates
where base_ave1<0.1 group by genotype, depth; quit;
proc sql; select distinct genotype, depth, mean(pred) as mean_pred from estimates
where base_ave1>=0.1 group by genotype, depth; quit;
proc sql; select distinct genotype, depth, mean(pred) as mean_pred, mean(StdErrPred) as MeanSTD_pred from estimates
 group by genotype, depth; quit;


proc glimmix data=dopa plots=all;
class animal genotype (ref='0') depth (ref='0');
model response=genotype|depth depth|base_ave1/solution or;
random animal/type=vc;
estimate 'Genotype' genotype 1 -1/exp; 
estimate 'Depth' depth 1 -1/exp; 

run;




