options formdlim=' ' nodate nocenter;
title;

* include macros to run algorithm;
%include '\\sfa18.partners.org\gogr$\Aim 1\R01_starting2013Mar\algorithm packages\SAS\ver3\crn-vdw\rec.vdw.macros.sas';
libname outdir '\\sfa18.partners.org\gogr$\Aim 1\R01_starting2013Mar\algorithm packages\SAS\ver3\crn-vdw\';

* import example datasets;
proc import out=exampledata 
    datafile= "X:\Aim 1\R01_starting2013Mar\algorithm packages\SAS\ver1\crn-vdw\exampledata_crnvdw.xls" 
	dbms=xls replace;
	getnames=yes;
run;

* lung cancer algorithm using default cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg);
data outdir.lungoutput;
 set combinedresults;
run;

* lung cancer algorithm using specified cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='lung', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.20);

* lung cancer algorithm using default cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg);

* lung cancer algorithm using specified cutoff to classify patients as recurrent;
%recvdw(inputdata=exampledata, disease='crc', ptid=case, month=month,
				   sec=nsmnoln, chm=nchm, hospice=nhsp, image=nimg, cutoff=0.10);
data outdir.crcoutput;
 set combinedresults;
run;

