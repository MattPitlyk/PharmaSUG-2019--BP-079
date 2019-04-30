/*******************************************************************************
| Program Name   : sas_fuzzy_match.sas
| Program Version: 1.0
| Program Purpose: Supports Fuzzy Matching algorithms in SAS 
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
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "&path.\copd_&src_data..csv"
       delimiter = '^' 
       MISSOVER DSD lrecl=32767 firstobs=1;

   informat cmindc $500. ;
   informat copdfl $500. ;
   format   cmindc $500. ;
   format   copdfl $500. ;
   input    cmindc $     
            copdfl $     ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
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

%macro FuzzyMatch(datain=&src_data._unique
                 ,dataout=FUZZY
                 ,word_list=
                 ,percent=
                 ,scoretype=
                 ,termnum=
                 );

*** Count number of words in proxy list ***;
data dum1;
   varx  = strip("&word_list.");
   varx_ = compress(varx);
   diffx = length(varx) - length(varx_) + 1;
run;

%global num_words;
proc sql noprint;
   select max(diffx) into : num_words from dum1;  *** Number of words in WORD_LIST ***;
quit;


data ds1;
   set &datain.;

   SUBSTRING = "&word_list.";

*** Remove separators in CMINDC which may create separate words ***;
   _CMINDX  = tranwrd(cmindc,"/"," ");
   _CMINDX  = tranwrd(_CMINDX,"-"," ");
   _CMINDX  = tranwrd(_CMINDX,";"," ");
   _CMINDX  = tranwrd(_CMINDX,","," ");
   CMINDX  = compbl(_CMINDX);
   CMINDX_ = compress(CMINDX);
   diffx = length(CMINDX) - length(CMINDX_) + 1; *** Number of words in CMINDX ***;
   nsubstr = diffx + 1 - &num_words.;            *** Number of substrings to check ***;
   drop CMINDX_ diffx _CMINDX;
run;

%global nsubstr;
proc sql noprint;
   select max(nsubstr) into : nsubstr from ds1;  *** Max Number of substrings to check ***;
quit;

*** Parse CMINDC into separate Substrings ***;
data ds2;
   set ds1;
   starttemp=cmindx;
   TEMP_WL = CMINDX;
   length CMINDC1-CMINDC%eval(&nsubstr.) $100;
   array trmlst{%eval(&nsubstr.)} $ CMINDC1-CMINDC%eval(&nsubstr.);
   do ii = 1 to nsubstr; 
      endpos = 1;
      do jj = 1 to &num_words.;
         if jj = 1 then do;
            endpos = find(TEMP_WL," ");
            newpos = endpos;
         end;
         else do;
            endpos = find(TEMP_WL," ",newpos+1);
            newpos = endpos;
         end;
      end;
      trmlst{ii} = strip(substr(TEMP_WL,1,endpos));
      if ii < nsubstr then TEMP_WL = strip(substr(TEMP_WL,find(TEMP_WL," ")+1));
   end;
   drop ii temp_wl jj newpos endpos;
run;


*** Pairwise, compare each proxy term with each word in CMINDC ***;
data ds3; 
   set ds2;

   if SUBSTRING ne "" then SOUNDEX0 = SOUNDEX(SUBSTRING);
   else SOUNDEX0 = "";

   array trmlst{%eval(&nsubstr.)} $ CMINDC1-CMINDC%eval(&nsubstr.);
   array sndxcm{%eval(&nsubstr.)} $ SNDXCM1-SNDXCM%eval(&nsubstr.);
   do kk = 1 to nsubstr;
      if trmlst{kk} ne "" then sndxcm{kk} = SOUNDEX(trmlst{kk});
      else sndxcm{kk} = "";
   end;

   array scores{%eval(&nsubstr.)}   SCORE1-SCORE%eval(&nsubstr.);
   do ii = 1 to nsubstr;
      if cmiss(SUBSTRING,trmlst{ii}) = 0 then do;
         _SCORE1 = COMPGED(SUBSTRING,trmlst{ii});         /* ScoreType 1 */
         _SCORE2 = SPEDIS(SUBSTRING ,trmlst{ii});         /* ScoreType 2 */
         _SCORE3 = COMPGED(SOUNDEX0 ,sndxcm{ii});         /* ScoreType 3 */
         _SCORE4 = SPEDIS(SOUNDEX0  ,sndxcm{ii});         /* ScoreType 4 */
         _SCORE5 = mean(_SCORE1,_SCORE2,_SCORE3,_SCORE4); /* ScoreType 5 */
         scores{ii} = _SCORE&scoretype.; 
      end;
      else scores{ii} = .;
   end;
 
   drop ii kk;
run;

*** Calculate User specified percent ***;
data ds3_temp(keep=score);
   set ds3(keep=SCORE: nsubstr);
   array scores{%eval(&nsubstr.)}   SCORE1-SCORE%eval(&nsubstr.);
   do ii = 1 to nsubstr;
      if scores{ii} ne . then do;
         SCORE = scores{ii};
         output;
      end;
   end;
run;   

proc sql noprint;
   select min(score) + &percent.*(max(score)-min(score)) into : pcntl from ds3_temp; 
quit;

data ds4; 
   set ds3;

   PERCENT = &percent.;
   array scores{%eval(&nsubstr.)}   SCORE1-SCORE%eval(&nsubstr.);

   POSSIBLE = 0;
   TRUEFL = 0;
   do ii = 1 to nsubstr;
      if . < scores{ii} <= &pcntl. then do;
         POSSIBLE+1;
         TRUEFL=1;
         end;
   end;

   %CalcStats(var=TRUEFL);

   keep CMINDC COPDFL COPDYN SUBSTRING CMINDX POSSIBLE PERCENT TRUEFL
        accuracy TYPEI_ERROR TYPEII_ERROR Precision Recall;
run;

proc sort data=ds4 out=&dataout.&termnum.;
   by descending POSSIBLE CMINDC;
run;
%mend FuzzyMatch;


proc format lib=work;
   value method 1 = "COMPGED"
                2 = "SPEDIS"
                3 = "COMPGED(SOUNDEX)"
                4 = "SPEDIS(SOUNDEX)"
                5 = "AVERAGE"
                ;
run;
quit;


%macro ChkStats(scoretype=,outdata=,ii=);

data labels;
   scoretype=&scoretype.;
   method = strip(put(scoretype,method.));
   THRESHOLD = &ii./100;
run;


   %put Checking Threshold = %eval(&ii.)%;

   options nonotes nosource;
   %FuzzyMatch(termnum=1,percent=&ii./100,scoretype=&scoretype.,word_list=COPD);
   %FuzzyMatch(termnum=2,percent=&ii./100,scoretype=&scoretype.,word_list=CHRONIC BRONCHITIS);
   %FuzzyMatch(termnum=3,percent=&ii./100,scoretype=&scoretype.,word_list=EMPHYSEMA);
   %FuzzyMatch(termnum=4,percent=&ii./100,scoretype=&scoretype.,word_list=CHRONIC OBSTRUCTIVE PULMONARY DISEASE);

   proc sort data=FUZZY1 out=FUZZY1a(keep=CMINDC COPDFL COPDYN TRUEFL); by CMINDC COPDFL COPDYN; run;
   proc sort data=FUZZY2 out=FUZZY2a(keep=CMINDC COPDFL COPDYN TRUEFL); by CMINDC COPDFL COPDYN; run;
   proc sort data=FUZZY3 out=FUZZY3a(keep=CMINDC COPDFL COPDYN TRUEFL); by CMINDC COPDFL COPDYN; run;
   proc sort data=FUZZY4 out=FUZZY4a(keep=CMINDC COPDFL COPDYN TRUEFL); by CMINDC COPDFL COPDYN; run;

   data FUZZY_ALL; 
      merge FUZZY1a(rename=(TRUEFL=T1))
            FUZZY2a(rename=(TRUEFL=T2)) 
            FUZZY3a(rename=(TRUEFL=T3)) 
            FUZZY4a(rename=(TRUEFL=T4));
      if T1=1 or T2=1 or T3=1 or T4=1 then TRUEFL=1;
      else TRUEFL=0;

      %CalcStats(var=TRUEFL);
   run;

   proc sql noprint;
      create table out2 as
      select &ii./100 format=6.2 as THRESHOLD , 
             mean(TRUEFL) as AVG_YES, mean(TYPEI_ERROR) as TYPEI_ERROR, mean(TYPEII_ERROR) as TYPEII_ERROR,
             mean(accuracy) as ACCURACY, mean(Precision) as PRECISION, mean(Recall) as RECALL, 
             sum(TYPEI_ERROR) as NumFalsePos, sum(TYPEII_ERROR) as NumFalseNeg
      from FUZZY_ALL;

   quit;
   run;

   data labels;
      merge labels
            out2; 
      by threshold;
   run;
   options notes source;

proc sort data=labels out=&outdata._&scoretype.;
   by threshold;
run;

%mend ChkStats;

%ChkStats(scoretype=1,outdata=COMPGED,ii=26);
%ChkStats(scoretype=3,outdata=COMPGEDX,ii=22);


******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
******************************************************;
