/*******************************************************************************
| Program Name   : train_v_test.sas
| Program Version: 1.0
| Program Purpose: Compare results on Training vs Test datasets
|                  for PharmaSUG 2019 paper BP-079
| SAS Version    : 9.4
| Created By     : Michael S. Rimler, GlaxoSmithKline
| Date           : 30Apr2019
| Input Data     : Training and Test Datasets (copd_train.csv, copd_test.csv)
********************************************************************************/
options nofmterr;

%let path=C:;

/* Read in Training data */
data work._train;
%let _EFIERR_ = 0; /* set the ERR0R detection macro variable */
infile "&path.\copd_train.csv"
       delimiter = '^' 
       MISSOVER DSD lrecl=32767 firstobs=1;

   informat cmindc $500. ;
   informat copdfl $500. ;
   format   cmindc $500. ;
   format   copdfl $500. ;
   input    cmindc $     
            copdfl $     ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERR0R detection macro variable */
run;

data train;
   set _train;
   copdyn = input(copdfl,best.);
run;

proc sort data=train out=train_unique nodupkey; 
   by copdyn cmindc; 
run;

proc freq data=train_unique; 
   tables copdfl; 
run;


/* Read in Test data */
data work._test;
%let _EFIERR_ = 0; /* set the ERR0R detection macro variable */
infile "&path.\copd_test.csv"
       delimiter = '^' 
       MISSOVER DSD lrecl=32767 firstobs=1;

   informat cmindc $500. ;
   informat copdfl $500. ;
   format   cmindc $500. ;
   format   copdfl $500. ;
   input    cmindc $     
            copdfl $     ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERR0R detection macro variable */
run;

data test;
   set _test;
   copdyn = input(copdfl,best.);
run;

proc sort data=test out=test_unique nodupkey; 
   by copdyn cmindc; 
run;

proc freq data=test_unique; 
   tables copdfl; 
run;



*** Merge Training and Test data (COPD only) ***;
proc sort data=train_unique out=train_unique_copd(keep=cmindc); 
   by cmindc; 
   where copdyn=1; 
run;

proc sort data=test_unique out=test_unique_copd(keep=cmindc);  
   by cmindc; 
   where copdyn=1; 
run;

data train_v_test;
   merge train_unique_copd(in=a)
         test_unique_copd (in=b);
   by cmindc;
   if a then trainfl = "Y";
   else trainfl = "N";

   if b then testfl = "Y";
   else testfl = "N";
run;

proc freq data=train_v_test; 
   tables trainfl*testfl; 
run;

/* *** Look at new terms in test ***; */
/*proc sort data=test_method2; by cmindc;  run;*/
/*data test_test; merge test_method2(in=a)*/
/*                      train_v_test(in=b);*/
/*                      by cmindc;*/
/*                      if trainfl="N" and testfl="Y";*/
/**/
/*run;*/


******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
