libname worklib "H:\SAS Course";

%let md = (year(date)*12 + month(date) - 1925*12 - 11);
%let yd = (year(date));
data y_ret (keep = date yd md mktrf rf d_p e_p svar b_m ntis eqis lty tms dfy infl ik cay);
		    retain date yd md mktrf rf d_p e_p svar b_m ntis eqis lty tms dfy infl ik cay; *reorder;
set worklib.predictors ;
	d_p = log(D12/Index);
	e_p = log(E12/Index);
	tms = lty - rf;
	dfy = BAA - AAA;
	md = &md;
	yd = &yd;
run;

data y_ret_combine (keep = date yd mktrf);
	set y_ret; 
run;

%macro run_OOS_combine(var, start, end);

%do T = &start %to &end;

	data y_ret_oos (keep = date yd mktrf &var);
		set y_ret; 
	run;

	data y_ret_oos2; set y_ret_oos; 
		if yd >= &T then mktrf = .;
	run;

	proc arima data = y_ret_oos2;
		identify var = mktrf crosscorr = &var nlag = 2 noprint;
		estimate input = (1$&var) method = ML noprint;
		Forecast lead = 1 id = yd out = out1 noprint;
	quit; 

	data out1 (keep = yd forecast); set out1; 
		if yd = &T then output;
	run;

	%if &T = &start %then
		%do;
			data y_ret_pred; set out1;
			run;
		%end;
	%else
		%do;
			proc append base = y_ret_pred  data = out1;
			quit; 
		%end;

%end;

proc sql; 
	create table temp as
 	select a.*, b.forecast as &var
	  from y_ret_combine as a left join y_ret_pred as b
      on a.yd = b.yd;
quit; 

data y_ret_combine;
	set temp;
run;

%mend run_OOS_combine;

Data List;
Input var$;
cards;
	rf	
	d_p
	e_p
	b_m
	svar
	ntis
	eqis
	lty
	tms
	dfy
	infl	
	ik
	cay 
	;
run;
proc printto log="desktop\temp.log";
run;

data _null_;
set List;
call execute('%run_OOS_combine(var=' ||var||  ',start = 1960, end = 2020)' );
run;

proc printto;
run;

data Eval_all;
run;
%macro evaluation(signal);
	%let k = 1;
	proc expand data=y_ret_combine out=y_ret_oos3;
			convert mktrf = forecast_0 / transformout=(CUAVE); 
	quit;

	data y_ret_oos3;
		set y_ret_oos3;
		forecast_1=&signal;
	run;

	data y_ret_oos3;
		set y_ret_oos3;
		forecast_0=lag(forecast_0);
		if forecast_1 = . then forecast_0 = .;
	run;

	data y_ret_oos3 (keep = date yd mktrf forecast_1 resid_1--resid_0_p);
	set y_ret_oos3; 
		resid_1 = mktrf - forecast_1;
		resid_0 = mktrf - forecast_0;
		resid_1_abs =  abs(resid_1);
		resid_0_abs =  abs(resid_0);
		resid_1_p   =  abs(resid_1/mktrf);
		resid_0_p   =  abs(resid_0/mktrf);
	if resid_1 ne . then output;
	run;

	proc sql; 
		create table Eval as
	 	select  count(resid_1) as T,
				avg(resid_1_abs)        as MAD_1,   avg(resid_0_abs)   as MAD_0,
				avg(resid_1_p)          as MAPE_1,  avg(resid_0_p)     as MAPE_0,
				avg(resid_1**2)         as MSE_1,   avg(resid_0**2)    as MSE_0,
				sqrt(CALCULATED MSE_1)  as RMSE_1,  sqrt(CALCULATED MSE_0)  as RMSE_0,
				(1 - (CALCULATED MSE_1/CALCULATED MSE_0)) as Rsq
	 	from y_ret_oos3;
	quit; 

	proc sql; 
		create table Eval2 as
	 	select MAD_0 - MAD_1 as MAD_diff,T,
			   MAPE_0 - MAPE_1 as MAPE_diff,
			   MSE_0  - MSE_1  as MSE_diff,
			   RMSE_0 - RMSE_1 as RMSE_diff,
			   Rsq, Rsq - (1-Rsq)*((&k )/(T-&k-1)) as adj_Rsq
	 	from Eval;
	quit;
	
	data Eval2;
		length var $ 10;
		set Eval2;
		var="&signal";
	run;

	data Eval_all;
		set Eval_all Eval2;
	run;
	%mend evaluation;

proc printto log="desktop\temp.log";
run;

data _null_;
set List;
call execute('%evaluation(' ||var||  ')' );
run;

proc printto;
run;

data Eval_all;
	set Eval_all;
	where var is not null;
run;

proc print data=Eval_all;
	title "Forecast for all signal";
run;

%macro SSE(start,end);
	data d;
		set y_ret_combine;
	run;
	data SSE;
	run;
	%do t= &start %to &end;
		data y_ret_combine;
			set d;
		run;

		data y_ret_combine;
			set y_ret_combine;
			where yd<=&t;

		data Eval_all;
		run;

		data _null_;
			set List;
			call execute('%evaluation(' ||var||  ')' );
		run;

		data Eval_all;
			set Eval_all;
			where var is not null;
		run;

		data Eval_all;
			set Eval_all;
			SSE_diff=MSE_diff*T;
		run;

		proc transpose data=Eval_all out=Eval_all;
			id var;
			var SSE_diff;
		run;

		data Eval_all;
			set Eval_all;
			drop _name_;
			end_year=&t;
		run;

		data SSE;
			set SSE Eval_all;
		run;

	%end;

	data SSE;
		set SSE;
		where end_year is not null;
	run;

	data y_ret_combine;
		set d;
	run;
%mend SSE;
proc printto log="desktop\temp.log";
run;

%SSE(1960,2020);

proc printto;
run;

proc sgplot data=sse;
	series x=end_year y=ik;
	title"cumulative SSE_diff for ik";
run;

%let signal=mean(of rf--cay);

proc expand data=y_ret_combine out=y_ret_oos3;
			convert mktrf = forecast_0 / transformout=(CUAVE); 
quit;

data y_ret_oos3;
	set y_ret_oos3;
	forecast_1=&signal;
run;

data y_ret_oos3;
	set y_ret_oos3;
	forecast_0=lag(forecast_0);
	if forecast_1 = . then forecast_0 = .;
run;

data y_ret_oos3 (keep = date yd mktrf forecast_1 resid_1--resid_0_p);
set y_ret_oos3; 
	resid_1 = mktrf - forecast_1;
	resid_0 = mktrf - forecast_0;
	resid_1_abs =  abs(resid_1);
	resid_0_abs =  abs(resid_0);
	resid_1_p   =  abs(resid_1/mktrf);
	resid_0_p   =  abs(resid_0/mktrf);
if resid_1 ne . then output;
run;

proc sql; 
	create table Eval as
 	select  count(resid_1) as T,
			avg(resid_1_abs)        as MAD_1,   avg(resid_0_abs)   as MAD_0,
			avg(resid_1_p)          as MAPE_1,  avg(resid_0_p)     as MAPE_0,
			avg(resid_1**2)         as MSE_1,   avg(resid_0**2)    as MSE_0,
			sqrt(CALCULATED MSE_1)  as RMSE_1,  sqrt(CALCULATED MSE_0)  as RMSE_0,
			(1 - (CALCULATED MSE_1/CALCULATED MSE_0)) as Rsq
 	from y_ret_oos3;
quit; 

%let k = 1;
proc sql; 
	create table Eval2 as
 	select MAD_0 - MAD_1 as MAD_diff,
		   MAPE_0 - MAPE_1 as MAPE_diff,
		   MSE_0  - MSE_1  as MSE_diff,
		   RMSE_0 - RMSE_1 as RMSE_diff,
		   Rsq, Rsq - (1-Rsq)*((&k )/(T-&k-1)) as adj_Rsq
 	from Eval;
quit;

proc print data=Eval2;
	title "Conbined forecast";
run;

proc sql;
CREATE TABLE fomc AS
SELECT
	a.*
	,b.a_flag AS flag
FROM
	worklib.mktrf_daily AS a left join worklib.fomc_dates AS b
on
	a.date=b.date;

data fomc;
	set fomc;
	retain cycle_a;
	if flag=1 then cycle_a=0;
	else cycle_a=cycle_a+1;
run;

proc sort data=fomc; 
	by descending date;
run;

data fomc;
	set fomc;
	retain cycle_d;
	if flag=1 then cycle_d=0;
	else cycle_d=cycle_d-1;
run;

data fomc;
	set fomc;
	if -6<=cycle_d<=-1 then cycle=cycle_d;
	else cycle=cycle_a;
run;

proc sort data=fomc;
	by date;
run;

data fomc;
	set fomc;
	where cycle is not null;
run;

data fomc;
	set fomc;
	if -1<=mod(cycle,10)<=3 or mod(cycle,10)=9 then w=1;
	else w=0;
run;


data fomc;
	set fomc;
	port=dm*w;
run;

proc reg data=fomc outest=reg noprint;
	model port=dm;
run;

proc sql;
CREATE TABLE stats AS
(SELECT
	"weighted" AS portfolio
	,MEAN(a.port)*252 AS ret
	,STD(a.port)*SQRT(252) AS std
	,CALCULATED ret/CALCULATED std AS sr
	,CALCULATED ret-0.5*4*CALCULATED std**2 AS cer
	,(MEAN(a.port-a.dm)*252)/(STD(a.port-a.dm)*SQRT(252)) AS ir
	,(b.Intercept*252)/(b._RMSE_*SQRT(252)) AS ar
FROM
	fomc AS a, reg AS b)
UNION
(SELECT
	"mkt" AS portfolio
	,MEAN(dm)*252 AS ret
	,STD(dm)*SQRT(252) AS std
	,CALCULATED ret/CALCULATED std AS sr
	,CALCULATED ret-0.5*4*CALCULATED std**2 AS cer
FROM
	fomc);

proc print data=stats;
	title "stats for US";
run;

data fomc;
	set fomc;
	port=dm_exUS*w;
run;

proc reg data=fomc outest=reg noprint;
	model port=dm_exUS;
run;


proc sql;
CREATE TABLE stats AS
(SELECT
	"weighted" AS portfolio
	,MEAN(a.port)*252 AS ret
	,STD(a.port)*SQRT(252) AS std
	,CALCULATED ret/CALCULATED std AS sr
	,CALCULATED ret-0.5*4*CALCULATED std**2 AS cer
	,(MEAN(a.port-a.dm_exUS)*252)/(STD(a.port-a.dm_exUS)*SQRT(252)) AS ir
	,(b.Intercept*252)/(b._RMSE_*SQRT(252)) AS ar
FROM
	fomc AS a, reg AS b)
UNION
(SELECT
	"mkt" AS portfolio
	,MEAN(dm_exUS)*252 AS ret
	,STD(dm_exUS)*SQRT(252) AS std
	,CALCULATED ret/CALCULATED std AS sr
	,CALCULATED ret-0.5*4*CALCULATED std**2 AS cer
FROM
	fomc);

proc print data=stats;
	title "stats for non-US";
run;

