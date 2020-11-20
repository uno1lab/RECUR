options formdlim=' ' nodate nocenter;
title;

* include macros to run algorithm;
%include 'recur-medicare.sas';

* import example datasets;
proc import out=exampledata 
    datafile= "exampledata.xls" 
	dbms=xls replace;
	getnames=yes;
run;

* lung cancer algorithm using default cutoff to classify patients as recurrent;
%recmedicare(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsm, rad=nrad, chm=nchm, hospice=nhsp, image=nimg);

* lung cancer algorithm using specified cutoff to classify patients as recurrent;
%recmedicare(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsm, rad=nrad, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.20);

* lung cancer algorithm using default cutoff to classify patients as recurrent;
%recmedicare(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsm, rad=nrad, chm=nchm, hospice=nhsp, image=nimg);

* lung cancer algorithm using specified cutoff to classify patients as recurrent;
%recmedicare(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsm, rad=nrad, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.10);


