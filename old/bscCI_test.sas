data testdat1;
	set ovad.testdata;
	y = n - x;
	p = x/n;
	z = probit(1-alpha/2);
run;

proc transpose data=testdat1 out=x_data(rename=(_NAME_=outcome COL1=count));
	by i;
	var x y;
run;

proc binomial data=x_data out=exact;
	by i;
	bi/bs;
	ou outcome;
	weight count;
run;

data exactci;
	set exact(keep=item value test by_val);

	if item in ('EST_PI' 'B_LCI_PI' 'B_UCI_PI');
	i=input(by_val, 8.);
run;

proc sort data=exactci;
	by test i;
run;

proc transpose data=exactci out=exactci;
	by test i;
	var value;
	id item;
run;

data final(keep=i Proportion LowerCL UpperCL method);
	set exactci;
	label Proportion = "Proportion"
	LowerCL = "95% Lower Confidence Limit"	
	UpperCL = "95% Upper Confidence Limit";
	Proportion = est_pi;
	LowerCL = b_lci_pi;
	UpperCL = b_uci_pi;
	length method $25.;
	method = "12. Blyth-Still-Casella";
run;
