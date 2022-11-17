/*** Code for cleaning/ creating new variables from Nigeria NAIIS Children's data ***
******* Dataset= NAIIS 2018 survey, children < 15 only********
******* Author: Colleen Leonard ***********
******* Date: 11/02/2021 ******************/

*Importing Excel file;
proc import datafile = "Nigeria_children_dataset_merged.xlsx"
	out = work.children
	dbms = xlsx
	replace
;
run; 

Proc contents data= children;
Run;

*Exploring the data;
Proc print data= children ;
	var age agegrp;
	where gst > 500;
RUN;
Proc print data= children ;
	var age ;
	where age= . OR age > 14;
RUN;
Proc print data= children ;
	var age ;
	where HRP2_1= .;
RUN;
*age is not missing for any child;
*agegrp is missing- need to re-create from Age;
*225 observations (0.7% of total) where gst > 500, cannot use their serology data;

Proc print data= children ;
	var PTID IND_ID HH_ID CHML CHMA age ;
	where age= 0 ;
RUN;

Data new_children;
	set children;

if gst >500 then delete; *225 samples, Comment out for the Antigen full dataset;

logAMA1 = log(ama1);
logCSP = log(csp);
logGLURP = log(glurp0);
logHRP2 = log(hrp2);
logLSA1 = log(lsa1);
logPfMSP1 = log(pfmsp1);
logPmMSP1 = log(pmmsp1);
logPoMSP1 = log(pomsp1);
logPvMSP1 = log(pvmsp1);

**Binary classification - of serology data ***;
if pfmsp1 gt 150 then PfMSP1pos=1; 
else PfMSP1pos=0; 

if AMA1 gt 209 then AMA1pos=1; 
else AMA1pos=0;
 
if GLURP0 gt 169 then GLURPpos=1; 
else GLURPpos=0; 

if CSP gt 362 then CSPpos=1; 
else CSPpos= 0; 

if HRP2 gt 337 then HRP2pos=1; 
else HRP2pos= 0; 

if LSA1 gt 124 then LSA1pos=1; 
else LSA1pos= 0; 

if pmmsp1 gt 550 then PmMSP1pos=1; 
else PmMSP1pos=0; 

if pomsp1 gt 362 then PoMSP1pos=1; 
else PoMSP1pos=0; 

if pvmsp1 gt 193 then PvMSP1pos=1; 
else PvMSP1pos=0; 

**Cutoffs for Antigen data ***;
if HRP2_1= . then HRP2_pos= .;
else if HRP2_1 gt 633 then HRP2_pos=1; 
else HRP2_pos=0; 

rename HRP2_1= HRP2_a;

if pAldo= . then pAldo_pos=.; 
else if pAldo gt 1865 then pAldo_pos=1; 
else pAldo_pos=0; 

if pLDH= . then pLDH_pos= .;
else if pLDH gt 430 then pLDH_pos=1; 
else pLDH_pos=0; 

if PvLDH= . then PvLDH_pos= .;
else if PvLDH gt 600 then PvLDH_pos=1; 
else PvLDH_pos=0; 

if age lt 5 and age ne . then agegrp= 1;
else if age lt 10 then agegrp= 2;
else if age lt 15 then agegrp= 3;

run;

*Exporting cleaned Ab dataset;
Proc export data= new_children (KEEP= ptid varstrat varunit A1QNumber HH_ID A1Q120 hhwgt AHMEMBER ab_wt ag_wt AHZone AHState AHLGA age agegrp gender AHTYPE wealthquintile A1Q119 A1Q115 A1Q108 A1Q109 A1Q105N 
								A1Q105A A1Q101A pfmsp1 pvmsp1 pomsp1 pmmsp1 ama1 glurp0 csp hrp2 lsa1 HRP2_a PfMSP1pos PmMSP1pos PoMSP1pos PvMSP1pos AMA1pos GLURPpos CSPpos HRP2pos LSA1pos HRP2_pos pLDH_pos PvLDH_pos pAldo_pos)
 		OUTFILE= "yourfilepath\Nigeria_children_sero_abbrev.xlsx" 
        DBMS=XLSX REPLACE;
     PUTNAMES=YES;
RUN;

*Creating dataset with Antigen data complete;
Data Antigen_data;
	set new_children;
	if PvLDH_pos= . OR pLDH_pos= . OR pAldo_pos= . OR HRP2_pos= . then delete;
RUN;
Proc freq data= Antigen_data ;
	table PvLDH_pos pLDH_pos pAldo_pos HRP2_pos; 
RUN;
*Exporting cleaned Antigen dataset;
Proc export data= Antigen_data (KEEP= ptid varstrat varunit A1QNumber HH_ID A1Q120 hhwgt AHMEMBER ab_wt ag_wt AHZone AHState AHLGA age agegrp gender AHTYPE wealthquintile A1Q119 A1Q115 A1Q108 A1Q109 A1Q105N 
								A1Q105A A1Q101A pfmsp1 pvmsp1 pomsp1 pmmsp1 ama1 glurp0 csp hrp2 lsa1 HRP2_a PfMSP1pos PmMSP1pos PoMSP1pos PvMSP1pos AMA1pos GLURPpos CSPpos HRP2pos LSA1pos HRP2_pos pLDH_pos PvLDH_pos pAldo_pos)
 		OUTFILE= "yourfilepath\Nigeria_children_Antigen_abbrev.xlsx" 
            DBMS=XLSX REPLACE;
     PUTNAMES=YES;
RUN;



*** FMM plots and outputs ***;
**for binary classification of IgG response;
ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logCSP = /k=2 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logGLURP = /k=2 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logHRP2 = /k=2 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logLSA1 = /k=2 dist=normal;
run;

*** non-Pf targets ***;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logPmMSP1 = /k=2 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logPoMSP1 = /k=2 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logPvMSP1 = /k=2 dist=normal;
run;


*** more than 2 components ****; 
ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logAMA1 = /k=6 dist=normal;
run;

ods graphics on;
proc fmm data=new_children plots=density (bins=100);;
model logPfMSP1 = /k=6 dist=normal;
run;

*** scatterplot of MFI-bg values ***; 

proc sgplot data=new_children;
	scatter x=pfmsp1 y=pvmsp1; 
run; 
proc sgplot data=new_children;
	scatter x=pfmsp1 y=pmmsp1; 
run; 

*Linear regression for MFI-bg values for msp1;
Proc reg data= new_children;
	model pvmsp1= pfmsp1; 
RUN;
*slope is effectively 0- no significant association in quantitative values;
*R2- very small= 0.002;

Proc reg data= new_children;
	model pmmsp1= pfmsp1; 
RUN;



