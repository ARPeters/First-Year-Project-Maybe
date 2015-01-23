PROC IMPORT OUT= WORK.DS 
            DATAFILE= "C:\Users\Psychometrics\Documents\dsAddicts.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc print;
run;

*Fitting Full Cox PH Model;
proc phreg data=ds;
model SurvivalTime*Status(0)= Clinic Prison Dose / ties=breslow;
OUTPUT out=phtest RESSCH= schClinic schPrison schDose;
run;

*Testing PH Assumption for each predictor;
data EventsOnly;
set phtest;
if Status=1;
run;

proc rank data=EventsOnly out=Ranked ties=mean;
var SurvivalTime;
ranks timerank;
run;

proc corr data=ranked NOSIMPLE;
VAR schClinic schPrison schDose;
with timerank;
run;

*Extended Cox Model;
*Heaviside function, f(x) =0 or 1 at interval r=90, event time = 428 weeks;
proc phreg data=ds;
model SurvivalTime*Status(0)= Clinic Clinic*HFR90 Prison Dose / ties=breslow;
IF SurvivalTime<428 THEN HFR90=0; ELSE HFR90=1;
run;

proc phreg data=ds;
model SurvivalTime*Status(0)= Clinic Clinic*HFR4080120 Prison Dose / ties=breslow;
IF SurvivalTime<168 THEN HFR4080120=0; ELSE IF SurvivalTIme<358 THEN HFR4080120=1; ELSE HFR4080120=2;
run;

