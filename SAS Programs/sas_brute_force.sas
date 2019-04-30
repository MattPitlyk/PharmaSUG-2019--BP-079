/*******************************************************************************
| Program Name   : sas_brute_force.sas
| Program Version: 1.0
| Program Purpose: Supports Brute Force SAS Exact String Matching 
|                  for classification of COPD conmed indications
|                  for PharmaSUG 2019 paper BP-079
| SAS Version    : 9.4
| Created By     : Michael S. Rimler, GlaxoSmithKline
| Date           : 30Apr2019
| Input Data     : Unique CMINDC dataset (training or test)
********************************************************************************/
options nofmterr;

%let path=C:;
%let src_data=train; *** test ***;

/* Read in data */
data work._&src_data.;
%let _EFIERR_ = 0; /* set the ERR0R detection macro variable */
infile "&path.\copd_&src_data..csv" 
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

data &src_data.;
   set _&src_data.;
   copdyn = input(copdfl,best.);
run;

proc sort data=&src_data. out=&src_data._unique nodupkey; 
   by copdyn cmindc; 
run;

%macro CalcStats(var=);
   *** Accuracy ***;
   if COPDYN=&var. then accuracy=1;
   else accuracy=0;

   *** Type I Error (false positives) ***;
   if COPDYN=0 and &var.=1 then TYPEI_ERROR = 1;
   else TYPEI_ERROR=0;

   *** Type II Error (false positives) ***;
   if COPDYN=1 and &var.=0 then TYPEII_ERROR = 1;
   else TYPEII_ERROR=0;

   *** Precision ***;
   if &var.=1 and COPDYN=1 then Precision = 1;
   else if &var.=1 then Precision=0;
   else Precision=.;

   *** Recall (percentage of true positives that got out of all positives) ***;
   if COPDYN=1 and &var.=1 then Recall = 1;
   else if COPDYN=1 then Recall=0;
   else Recall=.;
%mend CalcStats;

/* Brute Force Low-Tech */
data &src_data._method1;
   set &src_data._unique; 
   if find(cmindc,"COPD") + 
      find(cmindc,"CHRONIC BRONCHITIS") +
      find(cmindc,"EMPHYSEMA") +
      find(cmindc,"CHRONIC OBSTRUCTIVE PULMONARY DISEASE") > 0
      then COPDYN_METH1=1;
   else COPDYN_METH1=0;

   %CalcStats(var=COPDYN_METH1);
run;

proc freq data=&src_data._method1;
   tables accuracy TYPEI_ERROR TYPEII_ERROR Precision Recall;
run;

proc freq data=&src_data._method1;
   tables COPDYN*COPDYN_METH1;
run;



/* Brute Force 100% Accuracy */
data &src_data._method2(rename=(COPDYN_METH1=COPDYN_METH2));
   set &src_data._unique; 

   if find(cmindc,"COPD") +
      find(cmindc,"CHRONIC BRONCHITIS") +
      find(cmindc,"EMPHYSEMA") +
      find(cmindc,"CHRONIC OBSTRUCTIVE PULMONARY DISEASE") > 0
      then do;
         if find(compress(tranwrd(cmindc,"-"," ")),"NONCOPD")>0 then COPDYN_METH1=0;
         else if cmindc in ("ANXIETY DUE TO COPD"
                            "DESATURATION DURING EXCERCICE (COPD)"
                            "HYPERGLYCEMIA DURING COPD EXACERBATION"
                            "PANIC ATTACK SECONDARY TO COPD"
                            "SECONDARY INFECTION OF EMPHYSEMA"
                            "SPUTUM WITH CHRONIC OBSTRUCTIVE PULMONARY DISEASE") then COPDYN_METH1=0;
         else COPDYN_METH1=1;
      end;
   else if find(compress(tranwrd(cmindc,"-"," ")),"COPD") + 
           find(compress(tranwrd(cmindc,"."," ")),"COPD") >0 then COPDYN_METH1=1;
   else if min(find(cmindc,"CHRONIC"),1)   + min(find(cmindc,"OBSTRUCTIVE"),1) + 
           min(find(cmindc,"PULMONARY"),1) + min(find(cmindc,"DISEASE"),1) > 2 and
           find(cmindc,"HEART") = 0 then COPDYN_METH1=1;
   else if min(find(cmindc,"CHRONIC"),1)   + min(find(cmindc,"BRONHITIS"),1)  
                                           + min(find(cmindc,"BRONCHITIS"),1) >= 2 then COPDYN_METH1=1;
   else if length(CMINDC)<=4 and 
           min(find(cmindc,"C"),1) + min(find(cmindc,"O"),1) + 
           min(find(cmindc,"P"),1) + min(find(cmindc,"D"),1) > 2 and 
           cmindc not in ("APOC" "BPCO" "BPOC" "CORP" "DROP" "EPOC" "MPOC" "PBCO" "POCD") then COPDYN_METH1=1;
   else if find(cmindc,"SYMBICORT")>0 then COPDYN_METH1=1;
   else if find(cmindc,"AIRWAY OBSTRUCTION")>0 and 
      cmindc ne "UPPER AIRWAY OBSTRUCTION" then COPDYN_METH1=1;
   else if find(cmindc,"OPD")>0 then COPDYN_METH1=1;
   else if find(compress(tranwrd(cmindc,"-"," ")),"RESCUEMED")>0 and 
           find(cmindc,"ASTHMA") = 0 then COPDYN_METH1=1;
   else if cmindc in ("EMFYSEMA"
                      "EMPHYSIMIA"
                      "EMPHYSYMIA"
                      "CHOPN"
                      "COIPD"
                      "COPOD"
                      "COPPD"
                      "CHRINIC OBSTUCTIVE PULMONARY DISEASE"
                      "CHRONIC AIRWAYS DISEASE"
                      "CHRONIC OBSTRUCTIVE LUNG D/O"
                      "CHRONIC OBSTUCTIVE PULMONARY DISESE"
                      "COBD EXECERBATION"
                      "COPE ACUTE EXACERBATION"
                      "COPT EXECERBATION"
                      "EXACERBATION OF BRONCHITIS"
                      "RELIEVER INHALER"
                      "RESCUE DRUG FOR SHORTNESS OF BREATH"
                      "RESCUE INHALER"
                      ) then COPDYN_METH1=1;
   else if cmindc in ("ANXIETY AFTER HER HUSBANDS DEATH"
                      "HYPERTENSIÓN"
                      "RESCUE FOR"
                      ) then COPDYN_METH1=1;
   else COPDYN_METH1=0;
   %CalcStats(var=COPDYN_METH1);

run;

proc freq data=&src_data._method2;
   tables accuracy TYPEI_ERROR TYPEII_ERROR Precision Recall;
run;

proc freq data=&src_data._method2;
   tables COPDYN*COPDYN_METH2;
run;


******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
