/*
===========================================================================

Project:   RECUR

           A new algorithm for identifying cancer recurrence.

File:      recur-medicare.sas
Authors:   Hajime Uno, Angel Cronin, Nikki Carroll 
		   (email: huno@ds.dfci.harvard.edu)
Version:   1.00
Date:      22nd November 2020

Copyright: (C) 2020, Michael Hassett and Debra Ritzwoller's joint research 
           group at Dana-Farber Cancer Institute and Kaiser Permanente 
           Colorado.

           This software is free for non-commercial use. It may be copied,
           modified, and redistributed provided that this copyright notice
           is preserved on all copies. The intellectual property rights of
           the algorithms used reside with Dana-Farber Cancer Institute, 
           Michael Hassett's research group. 
           
           You may not use this software, in whole or in part, in support
           of any commercial product without the express consent of the
           author.

           There is no warranty or other guarantee of fitness of this
           software for any purpose. It is provided solely "as is".
           
           If you are pursuing commercial software development that 
           incorporates our software please contact Michael Hassett for
           information on licensing. 
           (email: Michael_Hassett@DFCI.HARVARD.EDU) 

===========================================================================
*/
*------start of: macro 1: getpeak--------------------------;
%macro getpeak(n, item, month, ptid);

	data peaks;
	 set phase2data (where=(&month.>0));
		by &ptid.;
		if first.&ptid. then item_sum=0;
			item_sum + &item.;
		item_avg= item_sum / &month.;
		item_avg_pre=lag(item_avg) + 1;
		if &month.=1 then item_avg_pre=1;
		if item_avg_pre > 0 then
			item_gap_pct=(&item. + 1 - item_avg_pre) / item_avg_pre;
 	keep &ptid. &item. &month. item_sum item_avg item_avg_pre item_gap_pct;
	run;

	proc sql;
		create table maxpct as
		select	in1.&ptid,
				max(item_gap_pct) as item_gap_pct_max,
				sum(&item.)>0 as anycode
		from peaks as in1
		group by in1.&ptid.;
		quit;
	data maxpct; 
	 set maxpct;
	 	where anycode>0;
	run;

	proc sql;
		create table peak&n. as
		select in1.&ptid, min(&month.) as peaktime&n.
		from peaks as in1, maxpct as in2
		where in1.&ptid.=in2.&ptid. and item_gap_pct_max>=0 and item_gap_pct=item_gap_pct_max
		group by in1.&ptid;
		quit;

%mend getpeak;
*--------end of: macro 1: getpeak--------------------------;


*------start of: macro 2: rec.medicare.lung---------------------;
%macro recmedicarelung(inputdata, ptid, month,
				   		sec, rad, chm, hospice, image, cutoff);

*--set paramater estimates;
	* phase 1, regression coefficients;
	%let b0=-3.47411573;	* intercept;
	%let b1=2.08231418;		* secondary malignancy w/ LN 2+;
	%let b2=1.95946431;		* radiation 2+;
	%let b3=3.88304468;		* hospice 1+;
	%let b4=0.04451343;		* imaging, continuous per year;

	* phase 2, -a-;
	%let a1=0.22;			* secondary malignancy w/ LN;
	%let a2=1.21;			* chemotherapy;
	%let a3=-0.94;			* imaging;
	* phase 2, weights;
	%let w1=0.279;			* secondary malignancy w/ LN;
	%let w2=0.296;			* chemotherapy;
	%let w3=0.425;			* imaging;
	
*--apply blackout windows for chemotherapy and radiation therapy;
	data afterblackout;
	 set &inputdata.;
	 	if &month.<=6 then do;
			&rad.=0;
			&chm.=0;
		end;
	run;

*--phase 1;
	proc sql;
		create table phase1results as
		select	in1.&ptid.,
				max(1,max(&month.)) as followup_months,
				sum(&sec.)>=2 as seccat,
				sum(&rad.)>=2 as radcat,
				sum(&hospice.)>=1 as hospicecat,
				sum(&image.) / calculated followup_months * 12 as imageyr,
				&b0. + 
				&b1. * calculated seccat +
				&b2. * calculated radcat +
				&b3. * calculated hospicecat +
				&b4. * calculated imageyr as score,
				exp(calculated score) / (1 + exp(calculated score)) as prob,
				calculated prob > &cutoff. as recurrence
		from afterblackout as in1
		group by in1.&ptid.;
		quit;

*--phase 2;
	data phase2data;
		merge	afterblackout
				phase1results (in=rec where=(recurrence=1));
			by &ptid.; 
			if rec=1;
	run;
	proc sort data=phase2data; by &ptid. &month.; run;

	%getpeak(n=1, item=&sec., month=&month., ptid=&ptid.);
	%getpeak(n=2, item=&chm., month=&month., ptid=&ptid.);
	%getpeak(n=3, item=&image., month=&month., ptid=&ptid.);

	data integrate;
	 merge peak1 peak2 peak3;
		by &ptid.;
		pred1 = peaktime1 - &a1.;
		pred2 = peaktime2 - &a2.;
		pred3 = peaktime3 - &a3.;
		if pred1 ne . and pred2 ne . and pred3 ne . 
			then predtime = (pred1*&w1. + pred2*&w2. + pred3*&w3.) /
							(&w1. + &w2. + &w3.);
		else if pred1 ne . and pred2 ne .
			then predtime = (pred1*&w1. + pred2*&w2.) /
							(&w1. + &w2.);
		else if pred1 ne . and pred3 ne .
			then predtime = (pred1*&w1. + pred3*&w3.) /
							(&w1. + &w3.);
		else if pred2 ne . and pred3 ne .
			then predtime = (pred2*&w2. + pred3*&w3.) /
							(&w2. + &w3.);
		else predtime=max(pred1, pred2, pred3);		
	run;
	data phase2results;
	 merge	phase1results (in=rec where=(recurrence=1))
			integrate;
		by &ptid.;
		if rec=1;
		if predtime=. then predtime = followup_months / 2;
		keep &ptid. followup_months peaktime1 peaktime2 peaktime3 pred1 pred2 pred3 predtime;
	run;

*--combine results from phase I and phase II;
	data combinedresults;
	 merge phase1results
		   phase2results;
		by &ptid.;
		keep &ptid. prob recurrence predtime;
	run;
	proc sort data=combinedresults; by &ptid.; run;

*--some clean up;
	proc datasets;
		delete afterblackout chmpeak imgpeak integrate maxpct peaks phase2data secpeak peak1-peak3 maxpct
			   phase1results phase2results;
	quit;
%mend recmedicarelung;
*--------end of: macro 2: rec.medicare.lung---------------------;

*------start of: macro 3: rec.medicare.crc-----------------------;
%macro recmedicarecrc(inputdata, ptid, month,
				   		sec, rad, chm, hospice, image, cutoff);

*--set paramater estimates;
	* phase 1, regression coefficients;
	%let b0=-4.5174690;		* intercept;
	%let b1=3.0040192;		* secondary malignancy w/ LN 2+;
	%let b2=1.5806104;		* radiation 2+;
	%let b3=3.2047289;		* hospice 1+;
	%let b4=0.1123489;		* imaging, continuous per year;

	* phase 2, -a-;
	%let a1=-0.07;			* secondary malignancy w/ LN;
	%let a2=3.69;			* chemotherapy;
	%let a3=-1.08;			* imaging;
	* phase 2, weights;
	%let w1=0.387;			* secondary malignancy w/ LN;
	%let w2=0.500;			* chemotherapy;
	%let w3=0.113;			* imaging;
	
*--apply blackout windows for chemotherapy and radiation therapy;
	data afterblackout;
	 set &inputdata.;
	 	if &month.<=12 then do;
			&rad.=0;
			&chm.=0;
		end;
	run;

*--phase 1;
	proc sql;
		create table phase1results as
		select	in1.&ptid.,
				max(1,max(&month.)) as followup_months,
				sum(&sec.)>=2 as seccat,
				sum(&rad.)>=2 as radcat,
				sum(&hospice.)>=1 as hospicecat,
				sum(&image.) / calculated followup_months * 12 as imageyr,
				&b0. + 
				&b1. * calculated seccat +
				&b2. * calculated radcat +
				&b3. * calculated hospicecat +
				&b4. * calculated imageyr as score,
				exp(calculated score) / (1 + exp(calculated score)) as prob,
				calculated prob > &cutoff. as recurrence
		from afterblackout as in1
		group by in1.&ptid.;
		quit;

*--phase 2;
	data phase2data;
		merge	afterblackout
				phase1results (in=rec where=(recurrence=1));
			by &ptid.; 
			if rec=1;
	run;
	proc sort data=phase2data; by &ptid. &month.; run;

	%getpeak(n=1, item=&sec., month=&month., ptid=&ptid.);
	%getpeak(n=2, item=&chm., month=&month., ptid=&ptid.);
	%getpeak(n=3, item=&image., month=&month., ptid=&ptid.);

	data integrate;
	 merge peak1 peak2 peak3;
		by &ptid.;
		pred1 = peaktime1 - &a1.;
		pred2 = peaktime2 - &a2.;
		pred3 = peaktime3 - &a3.;
		if pred1 ne . and pred2 ne . and pred3 ne . 
			then predtime = (pred1*&w1. + pred2*&w2. + pred3*&w3.) /
							(&w1. + &w2. + &w3.);
		else if pred1 ne . and pred2 ne .
			then predtime = (pred1*&w1. + pred2*&w2.) /
							(&w1. + &w2.);
		else if pred1 ne . and pred3 ne .
			then predtime = (pred1*&w1. + pred3*&w3.) /
							(&w1. + &w3.);
		else if pred2 ne . and pred3 ne .
			then predtime = (pred2*&w2. + pred3*&w3.) /
							(&w2. + &w3.);
		else predtime=max(pred1, pred2, pred3);		
	run;
	data phase2results;
	 merge	phase1results (in=rec where=(recurrence=1))
			integrate;
		by &ptid.;
		if rec=1;
		if predtime=. then predtime = followup_months / 2;
		keep &ptid. followup_months peaktime1 peaktime2 peaktime3 pred1 pred2 pred3 predtime;
	run;

*--combine results from phase I and phase II;
	data combinedresults;
	 merge phase1results
		   phase2results;
		by &ptid.;
		keep &ptid. prob recurrence predtime;
	run;
	proc sort data=combinedresults; by &ptid.; run;

*--some clean up;
	proc datasets;
		delete afterblackout chmpeak imgpeak integrate maxpct peaks phase2data secpeak peak1-peak3 maxpct
			   phase1results phase2results;
	quit;
%mend recmedicarecrc;
*------end of: macro 3: rec.medicare.crc-------------------------;


*------start of: macro 4: rec.medicare--------------------------;
%macro recmedicare(inputdata, disease, ptid, month,
				   sec, rad, chm, hospice, image, cutoff=.);
	%if &disease.='lung' %then %do;
		%if &cutoff.=. %then %do;
			%let cutoff=0.702;
		%end;
		%recmedicarelung(inputdata=&inputdata., ptid=&ptid., month=&month.,
				   sec=&sec., rad=&rad., chm=&chm., hospice=&hospice., image=&image., cutoff=&cutoff.);
	%end;

	%if &disease.='crc' %then %do;
		%if &cutoff.=. %then %do;
			%let cutoff=0.400;
		%end;
		%recmedicarecrc(inputdata=&inputdata., ptid=&ptid., month=&month.,
				   sec=&sec., rad=&rad., chm=&chm., hospice=&hospice., image=&image., cutoff=&cutoff.);
	%end;
%mend recmedicare;
*--------end of: macro 4: rec.medicare--------------------------;
