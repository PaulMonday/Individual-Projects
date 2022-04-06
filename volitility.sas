libname worklib "H:\SAS Course";
data ff;
	set worklib.ff_5factors_daily;
data vix;
	set worklib.vix_daily;
run;

%let ret=mktrf;

data ff;
	set ff;
	lret=log(1+&ret);
run;

proc sql;
CREATE TABLE rv AS
SELECT
	month(date) AS month
	,year(date) AS year
	,sum(lret) AS lret
	,sum(lret**2) AS rv_m
	,sqrt(12*CALCULATED rv_m) AS rv
FROM
	ff
GROUP BY
	CALCULATED year,CALCULATED month;

data vix;
	set vix(keep=date vix);
	year=year(date);
	month=month(date);
run;
proc sort data=vix;
	by descending date;
run;
proc sort data=vix nodupkey;
	by year month;
run;
data d;
	merge rv vix;
	by year month;
run;
data d;
	set d;
	where date is not null;
proc sgplot data=d;
	histogram rv;
	title "rv distribution,ret=&ret";
run;
proc sgplot data=d;
	histogram vix;
	title "vix distribution,ret=&ret";
run;

data d;
	set d;
	id=_N_;
run;
proc sql;
CREATE TABLE d AS
SELECT
	a.*
	,b.rv as ev
FROM
	d as a, d as b
WHERE
	a.id=b.id-1;
proc reg data=d;
	model ev=rv vix;
	title "linear regression y=ev, x=rv vix,ret=&ret";
	output out=d p=pred;
run;

proc sgplot data=d;
	histogram pred;
	title "volatility forecast distribution,ret=&ret";
run;
proc means data=d mean median;
	var pred;
	title "stats of forecast distribution,ret=&ret";
run;
data d0;
	set d;
run;

%macro get_ratio(v,r);

data d;
	set d0;
	w=&v/pred;
run;

proc sql;
CREATE TABLE d AS
SELECT
	b.lret
	,a.*
FROM
	d AS a, d AS b
WHERE
	a.id=b.id-1;
data d;
	set d;
	excr=lret*w;
	ann_ex=excr*12;
	ann_lret=lret*12;
	ar=ann_ex-ann_lret;
run;
proc reg data=d outest=reg noprint;
	model ann_ex=ann_lret;
run;
proc sql;
CREATE TABLE stats AS
SELECT
	a.intercept as alpha
	,a._RMSE_ AS RMSE
	,b.*
FROM
	reg AS a
	,(SELECT
		mean(ann_ex) AS ex_mean
		,std(excr)*sqrt(12) AS ex_std
		,mean(ar) AS ar_mean
		,std(excr-lret)*sqrt(12) AS ar_std
	FROM
		d) AS b;

data stats;
	set stats;
	sharpe_ratio=ex_mean/ex_std;
	information_ratio=ar_mean/ar_std;
	utility=ex_mean-0.5*&r*ex_std**2;
	appraisal_ratio=alpha/RMSE;
run;
%mend;

%get_ratio(0.12,4);
proc sgplot data=d;
	histogram w;
	title "distribution of market portfolio weight,ret=&ret";
run;
proc print data=stats;
	var ex_mean ex_std sharpe_ratio information_ratio utility appraisal_ratio;
	title "stats for Q7,ret=&ret";
	run;

data d;
	set d0;
	w=1;
run;
proc sql;
CREATE TABLE d AS
SELECT
	b.lret
	,a.*
FROM
	d AS a, d AS b
WHERE
	a.id=b.id-1;
data d;
	set d;
	excr=lret*w;
	ann_ex=excr*12;
	ann_lret=lret*12;
	ar=ann_ex-ann_lret;
run;
proc reg data=d outest=reg noprint;
	model ann_ex=ann_lret;
run;
proc sql;
CREATE TABLE stats AS
SELECT
	a.intercept as alpha
	,a._RMSE_ as RMSE
	,b.*
FROM
	reg AS a
	,(SELECT
		mean(ann_ex) AS ex_mean
		,std(excr)*sqrt(12) AS ex_std
		,mean(ar) AS ar_mean
		,std(excr-lret)*sqrt(12) AS ar_std
	FROM
		d) AS b;

data stats;
	set stats;
	sharpe_ratio=ex_mean/ex_std;
	information_ratio=ar_mean/ar_std;
	utility=ex_mean-0.5*4*ex_std**2;
	appraisal_ratio=alpha/RMSE;
run;
proc print data=stats;
	var ex_mean ex_std sharpe_ratio information_ratio utility appraisal_ratio;
	title "stats for w=1,ret=&ret";
	run;

%macro get_ratio_loop(start,end,r);
data stats_all;
run;
%do i= &start %to &end;
	%get_ratio(&i/100,&r)
	data stats;
		set stats;
		target=&i;
		r=&r;
	run;
	data stats_all;
		set stats_all stats;
	run;
%end;
data stats_all;
	set stats_all;
	where r is not null;
run;
%mend;
%get_ratio_loop(6,22,4);
proc sgplot data=stats_all;
	series y=ex_std x=target;
	series y=sharpe_ratio x=target;
	series y=information_ratio x=target;
	series y=utility x=target;
	series y=appraisal_ratio x=target;
	title "ratio for r=6% to 22%,r=4,ret=&ret";
run;

*Q9;
%get_ratio_loop(6,22,7);
proc sgplot data=stats_all;
	series y=ex_std x=target;
	series y=sharpe_ratio x=target;
	series y=information_ratio x=target;
	series y=utility x=target;
	series y=appraisal_ratio x=target;
	title "ratio for r=6% to 22%,r=7,ret=&ret";
run;

