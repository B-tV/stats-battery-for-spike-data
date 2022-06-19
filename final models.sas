*
combined baselines
	logbase using catGT*catDV2inv (nothing sig.)
	logbase for D1 cells (sig.)
l-DOPA
	polaritynum using catDV2inv*logbase (sig.)
	logpostprerat_abs using 3-way interaction (sig.)
	logpostprerat_abs on D1 cells (t-test)
D1ag
	polaritynum interaction (sig.)
	logpostprerat_abs interaction (sig.)
	logpostprerat_abs on D1 cells (t-test)
;

*	DATA for BASELINES;
PROC IMPORT OUT= LDOPAWOI.BASELINESWOINJ 
            DATAFILE= "C:\Users\Zhoulab\Documents\My SAS Files\baselines
4SAS.xls" 
            DBMS=EXCELCS REPLACE;
     RANGE="Sheet1$"; 
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;
data baselineswoinj; 
	set lDOPAwoi.baselineswoinj(firstobs=1 obs=198);
run;


* 	BASELINES;
* interaction makes not much difference to metrics (i.e. main effects model ~ same metrics);
*NO significant covtests;
proc GLIMMIX data=baselineswoinj plot=all;
class catGT catDV2inv animal;
model logbase = catGT|catDV2inv /
	solution ddfm=kr2;
output out = logbasepred pred std;
random _residual_ / subject=animal type=cs;
covtest DIAGR; 
covtest homogeneity;
lsmeans catGT*catDV2inv / diff=control('WT' 'dorsal') adjust=sidak stepdown adjdfe=row e cl;
run; *grouping 	BIC		sig. DF		sig. p	notes
VC				476.37				>>.05	
CS(~.07)		476.77				>>.05	covtest ~.1
VC by...									
catGT 			478.80				>>.05
catDV2inv 		477.39				>>.05  	covtest ~.15
CS by...									
catGT			481.40				>>.05
	.03 Pitx, 	.07 WT
catDV2inv  		481.05				>>.05
	.04 dorsal, .04 non-dorsal
;

* for communicating dataset;
proc tabulate data = baselineswoinj;
class catGT catDV2inv binombase;
var logbase;
table catGT*logbase, n mean stddev;
table catDV2inv*logbase, n mean stddev;
table catGT*catDV2inv*logbase, n mean stddev;
table catGT*catDV2inv*binombase*logbase, n mean stddev; *binombase to ~ binarize logbase;
run;
* for communicating model;
proc tabulate data = logbasepred;
class catGT catDV2inv binombase;
var pred stderr;
table catGT*(pred stderr), n mean stddev;
table catDV2inv*(pred stderr), n mean stddev;
table catGT*catDV2inv*(pred stderr), n mean stddev;
*table catGT*catDV2inv*binombase*(pred stderr), n mean stddev; *binombase to ~ binarize logbase;
run; * CONCLUSION = no difference in baselines between any factor or interaction;


*	DATA for D1cellbaselines;
PROC IMPORT OUT= LDOPAWOI.D1cellbaselines 
            DATAFILE= "C:\Users\Zhoulab\Documents\My SAS Files\D1cell baselines.xls" 
            DBMS=EXCELCS REPLACE;
     RANGE="Sheet1$"; 
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;
data D1cellbaselines; 
	set lDOPAwoi.D1cellbaselines;
	if base_ave1 >= .2 then do; binombasenum = 1; end;
	else if base_ave1 < .2 then do; binombasenum = 0; end;
run;

* testing association between logbase and catGT;
proc glimmix data=D1cellbaselines plot=all maxopt=100;*abspconv=1e-6;
class catGT binombase animal;
model logbase = catGT / solution ddfm=kr2;
random _residual_ / subject=animal type=vc group=catGT;
output out = GLIMbaseD1cells pred stderr;
covtest DIAGR;
covtest homogeneity;
LSMEANS catGT / diff e cl ;
run; *		BIC				p's			notes
VC			28.29			.0162		
CS			29.66			~.09		diagr ~.5
VC\catGT	24.97			.0115		homog ~.024
CS\catGT	28.03			~.09		homog ~.07, diagr ~.8
 .06 Pitx	.1 WT


* for communicating dataset;
proc tabulate data = D1cellbaselines;
class catGT binombase;
var logbase;
table catGT*logbase, n mean stddev;
table catGT*binombase*logbase, n mean stddev; *binombase to ~ binarize logbase;
run;
* for communicating results;
proc tabulate data = GLIMbaseD1cells;
class catGT binombase;
var pred stderr;
table catGT*(pred stderr), n mean stddev;
table catGT*binombase*(pred stderr), n mean stddev; *binombase to ~ binarize logbase;
run; *CONCLUSION = Pitx a.d. cells have ~a log unit lower baseline than WT
such that catGT is basically synonymous with logbase, so including it in later D1 models is redundant;




*	l-DOPA;
PROC IMPORT OUT= LDOPAWOI.LDOPAWOINJ 
            DATAFILE= "C:\Users\Zhoulab\Documents\My SAS Files\lDOPAwoin
jpartcrossedgrouped4SAS.xls" 
            DBMS=EXCELCS REPLACE;
     RANGE="Sheet1$"; 
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN; *SORTED DATA HERE WAS FIRST BY LOGBASE, THEN DEPTH, THEN ANIMAL, THEN GENOTYPE;
data ldopawoinj; 
	set lDOPAwoi.lDOPAwoinj(firstobs=1 obs=165);
	logpostprerat_abs=abs(logpostprerat);
	kreitzerdeltaave1_abs=abs(kreitzerdeltaave1);
	if logpostprerat >= 0 then do; polarity = 'up' ; polaritynum = 1; end;
	else if logpostprerat < 0 then do polarity = 'dn'; polaritynum= 0; end;
run;


*	POLARITYNUM
* having removed all insig. 2-ways... and before that having come from model with 3-way p.0563;
proc glimmix data=lDOPAwoinj plot=all maxopt=100 ; *where catGT='Pitx3Null' and logbase>-.79;
class catGT catDV2inv binombase polarity animal;
model polaritynum (ref=last) = catGT catDV2inv|logbase / 
	e3 oddsratio(unit logbase=-1 at logbase=-0.29 diff=all) solution ddfm=kr2 dist=binomial ;
random _residual_ / subject=animal type=cs group=catGT;
output out = lDOPApolaritylogbase pred stderr;
covtest DIAGR;
covtest homogeneity;

estimate
	'@dorsal logbase effect around mean' 
logbase [1, -1.29] [-1, -0.29]  
logbase*catDV2inv [1, -1.29 1] [-1, -0.29 1],
	'difference catDV2inv makes to logbase effect' 
catDV2inv*logbase [1, -1 1] [-1, -1 2]
	/ adjust=sidak stepdown adjdfe=row e CL exp ilink;

estimate 
		'logbase only' logbase [1, -1] [-1, 0],
		'catDV2inv AND base @ hi -0.29' catDV2inv*logbase [1, -0.29 1] [-1, -0.29 2],
		'catDV2inv AND base @ mean -0.79' catDV2inv*logbase [1, -0.79 1] [-1, -0.79 2],
		'catDV2inv AND base @ lo -1.29' catDV2inv*logbase [1, -1.29 1] [-1, -1.29 2],
		'catDV2inv with base @ lo-er -1.79' catDV2inv 1 -1 catDV2inv*logbase [1, -1.79 1] [-1, -1.79 2],
		'catDV2inv with base @ mean' catDV2inv 1 -1 catDV2inv*logbase [1, -0.79 1] [-1, -0.79 2],
		'catDV2inv with base @ hi -0.29' catDV2inv 1 -1 catDV2inv*logbase [1, -0.29 1] [-1, -0.29 2],
		'dorsal baseline difference (= O.R. estimate)' logbase [1, -1.29] [-1, -0.29] catDV2inv*logbase [1, -1.29 1] [-1, -0.29 1],
		'non-dorsal baseline difference(= -solution est)' logbase [1, -1.29] [-1, -0.29] catDV2inv*logbase [1, -1.29 2] [-1, -0.29 2],
		'difference baseline makes to depth effect (= row 4-row 2 = row 5-row 6 = row 8-row 9)' catDV2inv*logbase [1, -1 1] [-1, -1 2],
		'TYPE III catGT'	catGT 1 -1,
		'TYPE III catDV2inv' catDV2inv 1 -1
	/ e cl exp ilink;

LSMeans catDV2inv
	/ diff e at logbase=-1.29 cl odds or ilink;
LSMeans catDV2inv 
	/ diff e at logbase=-0.79 cl odds or ilink;
LSMeans catDV2inv 
	/ diff e at logbase=-0.29 cl odds or ilink;
LSMeans catDV2inv
	/ diff e bylevel cl odds or ilink;

run; * grouping	-2resLL	Gen.X2	GenX2/DF	p <.05				notes
CS(~.08)n/a		741.85	143.41	0.90		ALL but catDV2inv	catGT p ~.09, covtest ~.07		
catGT			735.27	160.00	1.00		logbase\*catDV2inv	diagr test ~.057
 .2 Pitx,	-.04 WT	
;

* for communicating dataset;
proc tabulate data = lDOPAwoinj;
class catGT catDV2inv binombase;
var logpostprerat_abs;
table catGT*logpostprerate_abs, n mean stddev;
table catDV2inv*logpostprerate_abs, n mean stddev;
table catGT*catDV2inv*logpostprerate_abs, n mean stddev;
table catGT*catDV2inv*binombase*logpostprerate_abs, n mean stddev; *binombase to ~ binarize logbase;
run;
* for communicating results;
proc tabulate data = lDOPApolaritylogbase;
class catGT catDV2inv binombase;
var pred stderr;
table catDV2inv*(pred stderr), n mean stddev;
table binombase*(pred stderr), n mean stddev;
table catDV2inv*binombase*(pred stderr), n mean stddev; 
*table catGT*catDV2inv*binombase*(pred stderr), n mean stddev; 
run; *CONCLUSION = logbase is significantly negatively correlated with polarity_num(=1???) with slope = -1.18 
(assuming mean effect of interaction [= .5*dorsal*logbase + .5*non-dorsal*logbase] 
added to logbase effect on its own, i.e. what's overall coded above as "logbase only"), 
AND with respect to depth, dorsal adds -1.64 units to the (already negative) logbase effect with respect to non-dorsal (the reference)
(i.e. the ln(OR) of logbase's [1-unit change] effect is decreased by 1.64 units in dorsal v. non-dorsal,
which ends up being a significant decrease, but also a significant overall group, 
i.e. dorsal*"low"base has significantly non-0 OR);


*  LOGPOSTPRERAT_ABS
* trying FULL model;
proc glimmix data = lDOPAwoinj plot=all; *where logbase < -0.5 ;
class catGT catDV2inv binombase polarity animal;
model logpostprerat_abs = catGT|catDV2inv|logbase / 
	e3 solution ddfm=kr2 ; 
random _residual_ / type=cs subject=animal group=catGT;
output out = lDOPAwoinjpred pred std;
covtest DIAGR;
covtest homogeneity;

estimate
	'catGT' catGT 1 -1,
	'logbase @-1.29' logbase -1.29,
	'logbase @-0.29' logbase -0.29,
	'logbase around mean' logbase [1, -1.29] [-1, -0.29],
	'@Pitx @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 1] [-1, -0.32 1] 
logbase*catDV2inv [1, -1.32 1] [-1, -0.32 1] 
catGT*catDV2inv*logbase [1, -1.32 1 1] [-1, -0.32 1 1],
	'@Pitx @non-dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 1] [-1, -0.32 1] 
logbase*catDV2inv [1, -1.32 2] [-1, -0.32 2] 
catGT*catDV2inv*logbase [1, -1.32 1 2] [-1, -0.32 1 2],
	'@WT @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 2] [-1, -0.32 2] 
logbase*catDV2inv [1, -1.32 1] [-1, -0.32 1] 
catGT*catDV2inv*logbase [1, -1.32 2 1] [-1, -0.32 2 1],
	'@WT @non-dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 2] [-1, -0.32 2] 
logbase*catDV2inv [1, -1.32 2] [-1, -0.32 2] 
catGT*catDV2inv*logbase [1, -1.32 2 2] [-1, -0.32 2 2],
	'difference catDV2inv makes to logbase effect @ WT'
logbase*catDV2inv [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 2 1] [-1, -1 2 2],
	'TYPE III catDV2inv' catDV2inv 1 -1,
	'TYPE III catGT*catDV2inv' catGT*catDV2inv 1 -1 -1 1
/ e cl;

estimate
	'@Pitx @dorsal logbase effect around mean' 
logbase [1, -1.29] [-1, -0.29] 
catGT*logbase [1, -1.29 1] [-1, -0.29 1] 
logbase*catDV2inv [1, -1.29 1] [-1, -0.29 1] 
catGT*catDV2inv*logbase [1, -1.29 1 1] [-1, -0.29 1 1],
	'difference catDV2inv makes to logbase effect @ Pitx'
logbase*catDV2inv [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 1 2],
	'difference catGT makes to logbase effect @ dorsal'
logbase*catGT [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 2 1]
	/ adjust=sidak stepdown adjdfe=row e CL;

LSMeans catGT*catDV2inv
	/ adjust=sidak stepdown adjdfe=row at logbase=-1.29 e cl;
LSMeans catGT*catDV2inv
	/ adjust=sidak stepdown adjdfe=row at logbase=-0.79 e cl;
LSMeans catGT*catDV2inv
	/ adjust=sidak stepdown adjdfe=row at logbase=-0.29 e cl;

LSMeans catGT
	/ diff at logbase=0 e cl;
LSMeans catGT
	/ diff bylevel e cl;

run; * grouping	BIC		DF			p <.05					notes
CS(~.08)n/a		281.12	21.52, 156.5catGT, logbase			hose-sprayed residuals
catGT			233.21	10.73, 143.4catGT, logbase, 3-way	
 .15 Pitx, .001 WT 
;

proc tabulate data = lDOPAwoinjpred;
class catGT catDV2inv binombase;
var pred stderr;
table catGT*(pred stderr), n mean stddev;
table binombase*(pred stderr), n mean stddev ;
table catGT*catDV2inv*binombase*(pred stderr), n mean stddev; 
run; *conclusion = ;


*	for dorsal D1 CELLS ONLY;
*filter dataset, firstobs-2, i.e. 167, will get the other two (non-dorsal) D1 cells;
data ldopawoinjD1cells; 
	set lDOPAwoi.lDOPAwoinj(firstobs=169 obs=180);
	logpostprerat_abs=abs(logpostprerat);
	kreitzerdeltaave1_abs=abs(kreitzerdeltaave1);
	if base_ave1 > 0.2 then binombase = '>0.2';
    else if base_ave1 <= 0.2 then binombase = '<=.2';
	if logpostprerat >= 0 then do; polarity = 'up' ; polaritynum = 1; end;
	else if logpostprerat < 0 then do polarity = 'dn'; polaritynum= 0; end;
run;


*ttest for D1cells according to mehmet;
proc ttest data=lDOPAwoinjD1cells cochran ci=equal umpu;
class catGT;
var logbase logpostprerat_abs kreitzerdeltaave1_abs;
run;




*	D1AG;
PROC IMPORT OUT= LDOPAWOI.D1agwoinj 
            DATAFILE= "C:\Users\Zhoulab\Documents\My SAS Files\D1agwoinj4SAS.xls" 
            DBMS=EXCELCS REPLACE;
     RANGE="Sheet1$"; 
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;
data D1agwoinj; 
	set lDOPAwoi.D1agwoinj(firstobs=1 obs=153);
	logpostprerat_abs=abs(logpostprerat);
	kreitzerdeltaave1_abs=abs(kreitzerdeltaave1);
	if logpostprerat >= 0 then do; polarity = 'up' ; polaritynum = 1; end;
	else if logpostprerat < 0 then do polarity = 'dn'; polaritynum= 0; end;
	if logbase <-0.82 then do; binomlogbase = '<-0.82' ; end;
	else if logbase >=-0.82 then do binomlogbase = '>-0.82' ; end;
run;


*	POLARITYNUM;
* full interactions;
proc glimmix data=D1agwoinj plot=all maxopt=100 ;*abspconv=1e-6;
*NLOPTIONS technique=nrridg ;*absfconv=1e-4 5 absgconv=1e-4; 
class catGT catDV2inv binombase polarity animal;
model polaritynum (ref='0') = catGT|catDV2inv|logbase / 
	e3 oddsratio(unit logbase=-1 at logbase=-0.32 diff=all) solution ddfm=kr2 dist=binomial ;
random _residual_ / subject=animal type=cs group=catDV2inv;

output out = D1agpolaritylogbase pred stderr;
covtest DIAGR;
covtest homogeneity;
estimate 
	'@Pitx @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 1] [-1, -0.32 1] 
logbase*catDV2inv [1, -1.32 1] [-1, -0.32 1] 
catGT*catDV2inv*logbase [1, -1.32 1 1] [-1, -0.32 1 1],
	'@Pitx @non-dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 1] [-1, -0.32 1] 
logbase*catDV2inv [1, -1.32 2] [-1, -0.32 2] 
catGT*catDV2inv*logbase [1, -1.32 1 2] [-1, -0.32 1 2],
	'@WT @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 2] [-1, -0.32 2] 
logbase*catDV2inv [1, -1.32 1] [-1, -0.32 1] 
catGT*catDV2inv*logbase [1, -1.32 2 1] [-1, -0.32 2 1],
	'@WT @non-dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 2] [-1, -0.32 2] 
logbase*catDV2inv [1, -1.32 2] [-1, -0.32 2] 
catGT*catDV2inv*logbase [1, -1.32 2 2] [-1, -0.32 2 2],
	'difference catDV2inv makes to logbase effect @ WT'
logbase*catDV2inv [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 2 1] [-1, -1 2 2],
	'TYPE III logbase' logbase -1,
	'TYPE III catGT' catGT 1 -1,
	'TYPE III catDV2inv' catDV2inv 1 -1,
	'TYPE III catGT*catDV2inv' catGT*catDV2inv 1 -1 -1 1,
	'difference catDV2inv makes @Pitx' 
catDV2inv 1 -1 
catGT*catDV2inv [1, 1 1] [-1, 1 2],
	'difference catDV2inv makes @WT' 
catDV2inv 1 -1
catGT*catDV2inv [1, 2 1] [-1, 2 2],	
	'difference catGT makes @dorsal' 
catGT 1 -1 
catGT*catDV2inv [1, 1 1] [-1, 2 1],
	'difference catGT makes @non-dorsal' 
catGT 1 -1
catGT*catDV2inv [1, 1 2] [-1, 2 2]
	/ e CL exp ilink;

estimate 
	'@Pitx @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*logbase [1, -1.32 1] [-1, -0.32 1] 
logbase*catDV2inv [1, -1.32 1] [-1, -0.32 1] 
catGT*catDV2inv*logbase [1, -1.32 1 1] [-1, -0.32 1 1],
	'difference catDV2inv makes to logbase effect @ Pitx'
logbase*catDV2inv [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 1 2],
	'difference catGT makes to logbase effect @ dorsal'
logbase*catGT [1, -1 1] [-1, -1 2] 
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 2 1]
	/ adjust=sidak stepdown adjdfe=row e CL exp ilink;

LSMeans 
	catGT*catDV2inv / e adjust=sidak stepdown adjdfe=row odds oddsratio CL ilink at logbase=-1.32;
LSMeans 
	catGT*catDV2inv / e adjust=sidak stepdown adjdfe=row odds oddsratio CL ilink at logbase=-0.82; 
LSMeans 
	catGT*catDV2inv / e adjust=sidak stepdown adjdfe=row odds oddsratio CL ilink at logbase=-0.32;
LSMeans 
	catGT*catDV2inv / e adjust=sidak stepdown adjdfe=row odds oddsratio CL ilink at logbase=-0.00;

LSMEstimate catGT*catDV2inv
	'catGT dorsal' [1, 1 1] [-1, 2 1],
	'Pitx catDV2inv' [1, 1 1] [-1, 1 2],
	'WT catDV2inv' [1, 2 1] [-1, 2 2]
	/ at logbase=-1.32 adjust=sidak stepdown adjdfe=row e CL exp ilink;
LSMEstimate catGT*catDV2inv
	'catGT dorsal' [1, 1 1] [-1, 2 1],
	'Pitx catDV2inv' [1, 1 1] [-1, 1 2],
	'WT catDV2inv' [1, 2 1] [-1, 2 2]
	/ at logbase=-0.82 adjust=sidak stepdown adjdfe=row e CL exp ilink;
LSMEstimate catGT*catDV2inv
	'catGT dorsal' [1, 1 1] [-1, 2 1],
	'Pitx catDV2inv' [1, 1 1] [-1, 1 2],
	'WT catDV2inv' [1, 2 1] [-1, 2 2]
	/ at logbase=-0.32 adjust=sidak stepdown adjdfe=row e CL exp ilink;
LSMEstimate catGT*catDV2inv
	'catGT dorsal' [1, 1 1] [-1, 2 1],
	'Pitx catDV2inv' [1, 1 1] [-1, 1 2],
	'WT catDV2inv' [1, 2 1] [-1, 2 2]
	/ at logbase=0.0 adjust=sidak stepdown adjdfe=row e CL exp ilink;

LSMEstimate catGT*catDV2inv
	'catGT dorsal' [1, 1 1] [-1, 2 1],
	'Pitx catDV2inv' [1, 1 1] [-1, 1 2]
	/ at logbase=0 adjust=sidak stepdown adjdfe=row e CL exp ilink;
LSMEstimate catGT*catDV2inv	
	'WT dorsal' [1, 2 1],
	'WT non-dorsal' [1, 2 2],
	'WT catDV2inv' [1, 2 1] [-1, 2 2]
	/ at logbase=0 e CL exp ilink;

run; * cov		-2resLL	gen.X2	X2/df	p<.05					notes
VC				686.20	148.24	1.02	catGT*catDV2inv\*logbasecatDV2inv p<.1
CS 																did not converge
catDV2inv		682.29	145.00	1.00	catGT*catDV2inv\*logbase
;

proc tabulate data = D1agpolaritylogbase;
class catGT catDV2inv binomlogbase;
var pred stderr;
table catGT*catDV2inv*(pred stderr), n mean stddev;
table catGT*catDV2inv*binomlogbase*(pred stderr), n mean stddev; 
run; *CONCLUSIONS: 
;


*	testing LOGPOSTPRERAT_ABS;
* trying mehmet's suggestion;
proc glimmix data = D1agwoinj plot=all maxopt=100;
class catGT catDV2inv binombase polarity animal;
model logpostprerat_abs = catGT catDV2inv logbase catGT*catDV2inv*logbase / 
	e3 solution ddfm=kr2; 
output out = D1agwoinjpred pred std;
random _residual_ / subject=animal type=cs group=catDV2inv;
covtest DIAGR; 
covtest homogeneity;

estimate 
	'logbase alone' 
logbase [1, -1.32] [-1, -0.32],
	'logbase @WT @dorsal' 
logbase [1, -1.32] [-1, -0.32] 
catGT*catDV2inv*logbase [1, -1.32 2 1] [-1, -0.32 2 1],
	'logbase @WT non-dorsal' 
logbase [1, -1.32] [-1, -0.32] 
catGT*catDV2inv*logbase [1, -1.32 2 2] [-1, -0.32 2 2],
	'difference catDV2inv makes to logbase effect @ WT'
catGT*catDV2inv*logbase [1, -1 2 1] [-1, -1 2 2],
	'logbase @Pitx non-dorsal' 
logbase [1, -1.32] [-1, -0.32] 
catGT*catDV2inv*logbase [1, -1.32 1 2] [-1, -0.32 1 2],
	'NOT Type III 3-way without 2-ways...' 
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 1 2] [-1, -1 2 1] [1, -1 2 2],
	'TYPE III catGT' catGT 1 -1,
 	'TYPE catDV2inv' catDV2inv 1 -1
	/ e CL;

estimate 
	'@Pitx @dorsal logbase effect around mean' 
logbase [1, -1.32] [-1, -0.32] 
catGT*catDV2inv*logbase [1, -1.32 1 1] [-1, -0.32 1 1],
	'difference catDV2inv makes to logbase effect @ Pitx'
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 1 2],
	'difference catGT makes to logbase effect @ dorsal'
catGT*catDV2inv*logbase [1, -1 1 1] [-1, -1 2 1]
	/ adjust=sidak stepdown adjdfe=row e CL;

LSMeans catDV2inv / at logbase=-1.32 e CL;
LSMeans catDV2inv / at logbase=-0.32 e CL;
LSMeans catGT / at logbase=-1.32 e CL;
LSMeans catGT / at logbase=-0.32 e CL;

run; *cov		BIC		DF for sig. p	p<.05			notes
CS by...
catDV2inv		237.62	127.6, 73.33	logbase\3-way	residuals a bit hosed down
 .1	dorsal	.01	non-dorsal
catGT*catDV2inv	230.45	112.2, 63.72	logbase\3-way	decent residuals
 .24	Pitx dorsal	.03	Pitx non-dorsal	.006	WT dorsal	.01	WT non-dorsal
;

proc tabulate data = D1agwoinjpred;
class catGT catDV2inv binomlogbase;
var pred stderr;
table catGT*catDV2inv*(pred stderr), n mean stddev;
table catGT*catDV2inv*binomlogbase*(pred stderr), n mean stddev; 
run; *CONCLUSIONS: 
;


*	for ONLY dorsal D1 cells (+2 to obs would give last two available, both non-dorsal);
*filter dataset (D1 cells from dorsal only since baselines seem to be different there, 
but no high baseline Pitx cells to test interaction with catGT);
data D1agwoinjD1cells; 
	set lDOPAwoi.D1agwoinj(firstobs=155 obs=165);
	logpostprerat_abs=abs(logpostprerat);
	kreitzerdeltaave1_abs=abs(kreitzerdeltaave1);
	if logpostprerat >= 0 then do; polarity = 'up' ; polaritynum = 1; end;
	else if logpostprerat < 0 then do polarity = 'dn'; polaritynum= 0; end;
run;


*ttest for D1cells according to mehmet;
proc ttest data=D1agwoinjD1cells cochran ci=equal umpu;
class catGT;
var logbase logpostprerat_abs kreitzerdeltaave1_abs;
run;


*	for plotting and/or saving separately...;
*for export must use replace option or delete xls's already in place...;
proc sql; 
create table ldopawoinjpred_saved 
as select distinct logpostprerat_abs, animal, catGT, catDV2inv, logbase, Pred, StdErr 
from ldopawoinjpred;
quit;
proc export 
data=ldopawoinjpred_saved 
outfile='C:\Users\Zhoulab\Documents\ben\misc figs-vids\my figures\ldopawoinjpred.saved_xlsx' dbms=xlsx;
run;
proc sql; 
create table ldopapolaritylogbase_saved 
as select distinct polaritynum, animal, catGT, catDV2inv, logbase, Pred, StdErr 
from ldopapolaritylogbase;
quit;
proc export 
data=ldopapolaritylogbase_saved 
outfile='C:\Users\Zhoulab\Documents\ben\misc figs-vids\my figures\ldopapolaritylogbase_saved.xlsx' dbms=xlsx;
run;

proc sql; 
create table D1agwoinjpred_saved
as select distinct logpostprerat_abs, animal, catGT, catDV2inv, logbase, Pred, StdErr 
from D1agwoinjpred;
quit;
proc export 
data=D1agwoinjpred_saved 
outfile='C:\Users\Zhoulab\Documents\ben\misc figs-vids\my figures\D1agwoinjpred_saved.xlsx' dbms=xlsx;
run;
proc sql;
create table D1agpolaritylogbase_saved
as select distinct polaritynum, animal, catGT, catDV2inv, logbase, Pred, StdErr 
from D1agwoinjpred;
quit;
proc export 
data=D1agpolaritylogbase_saved 
outfile='C:\Users\Zhoulab\Documents\ben\misc figs-vids\my figures\D1agpolaritylogbase_saved.xlsx' dbms=xlsx;
run;
