options formdlim=' ' nodate nocenter;
title;

* include macros to run algorithm;
%include 'recur-vdw.sas';
libname outdir '.';

* import example datasets;
proc import out=exampledata 
    datafile= "exampledata.xlsx" 
	dbms=xls replace;
	getnames=yes;
run;

*------------------------------------;
* Lung cancer algorithm using default cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg);

proc contents data=combinedresults ; run ;

*------------------------------------;
* Lung cancer algorithm using specified cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.20);

proc contents data=combinedresults ; run ;

*------------------------------------;
* Colorectal cancer algorithm using default cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg);

proc contents data=combinedresults ; run ;

*------------------------------------;
* Colorectal cancer algorithm using specified cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.10);

proc contents data=combinedresults ; run ;


