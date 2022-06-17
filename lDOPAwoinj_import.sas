PROC IMPORT OUT= LDOPAWOI.lDOPAwoinjpartcrossedgrouped 
            DATAFILE= "C:\Users\Zhoulab\Documents\My SAS Files\lDOPAwoin
jpartcrossedgrouped4SAS.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
