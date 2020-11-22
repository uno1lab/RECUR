/*
===========================================================================

Project:   RECUR

           A new algorithm for identifying cancer recurrence.

File:      recur-vdw-breast.sas
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


*------start of: macro 2: rec.vdw.breast-----------------------;
%macro recvdwbreast(inputdata, ptid, month,
				   		sec, chm, rad, hospice, mast, image, obs, cutoffph1);

*--set paramater estimates;
	* phase 1, regression coefficients;
	%let b0=-4.633634225;	* intercept;
	%let b1=2.58750348;		* secondary malignancy w/o LN 1;
	%let b2=4.874697311;	* secondary malignancy w/o LN 2;
	%let b3=4.966462159;	* secondary malignancy w/o LN 3+;
	%let b4=1.329849979;	* chemotherapy 2+;
	%let b5=1.661065818;	* radiation 2+;
	%let b6=0.819394667;	* hospice 1+;
	%let b7=0.882805651;	* mastectomy 2+;
	%let b8=0.1045132;		* imaging, continuous per year;

	* phase 2, -a-;
	%let a1=1.379;			* secondary malignancy w/o LN;
	%let a2=2.317;			* chemotherapy;
	%let a3=1.892;			* radiation;
	%let a4=0.501;			* mastectomy;
	%let a5=0.398;			* observation;
	* phase 2, weights;
	%let w1=0.213;			* secondary malignancy w/ LN;
	%let w2=0.273;			* chemotherapy;
	%let w3=0.297;			* imaging;
	%let w4=0.104;			* mastectomy;
	%let w5=0.113;			* observation;

	* phase 3, regression coefficients;
	%let ph3b0=-1.198822903;
	%let ph3b1=2.2839164;
	%let ph3b2=2.375099222;
	%let ph3b3=1.264097188;
	%let ph3b4=1.738575672;

*--apply blackout windows for treatment (chemotherapy, radiation, and mastectomy);
	data afterblackout;
	 set &inputdata.;
	 	if &month.<=12 then do;
			&chm.=0;
			&rad.=0;
			&mast.=0;
		end;
	run;

*--phase 1;
	proc sql;
		create table phase1calc as
		select	in1.&ptid.,
				max(1,max(&month.)) as followup_months,
				sum(&sec.) as totalsec,
				sum(&sec.)=1 as seccat1,
				sum(&sec.)=2 as seccat2,
				sum(&sec.)>=3 as seccat3,
				sum(&chm.)>=2 as chmcat,
				sum(&rad.)>=2 as radcat,
				sum(&hospice.)>=1 as hospicecat,
				sum(&mast.)>=2 as mastcat,
				sum(&image.) / calculated followup_months * 12 as imageyr,
				&b0. + 
				&b1. * calculated seccat1 +
				&b2. * calculated seccat2 +
				&b3. * calculated seccat3 +
				&b4. * calculated chmcat +
				&b5. * calculated radcat +
				&b6. * calculated hospicecat +
				&b7. * calculated mastcat +
				&b8. * calculated imageyr as score,
				exp(calculated score) / (1 + exp(calculated score)) as expscore
		from afterblackout as in1
		group by in1.&ptid.;
		quit;
	data phase1results;
	 set phase1calc;
	 	if totalsec > 34 then prob=1; else prob=expscore;
		if prob > &cutoffph1. then recurrence=1; else recurrence=0;
	run;

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
	%getpeak(n=3, item=&rad., month=&month., ptid=&ptid.);
	%getpeak(n=4, item=&mast., month=&month., ptid=&ptid.);
	%getpeak(n=5, item=&obs., month=&month., ptid=&ptid.);

	data integrate;
	 merge peak1 peak2 peak3 peak4 peak5;
		by &ptid.;

		pred1 = peaktime1 - &a1.;	product1 = pred1*&w1.;
		pred2 = peaktime2 - &a2.;	product2 = pred2*&w2.;
		pred3 = peaktime3 - &a3.;	product3 = pred3*&w3.;
		pred4 = peaktime4 - &a4.;	product4 = pred4*&w4.;
		pred5 = peaktime5 - &a5.;	product5 = pred5*&w5.;

		if peaktime1 ne . then weight1=&w1.;
		if peaktime2 ne . then weight2=&w2.;
		if peaktime3 ne . then weight3=&w3.;
		if peaktime4 ne . then weight4=&w4.;
		if peaktime5 ne . then weight5=&w5.;

		predtime = sum(product1, product2, product3, product4, product5) / 
				   sum(weight1, weight2, weight3, weight4, weight5);
	run;
	data phase2results;
	 merge	phase1results (in=rec where=(recurrence=1))
			integrate;
		by &ptid.;
		if rec=1;
		if predtime=. then predtime = followup_months / 2;
		keep &ptid. followup_months peaktime1 peaktime2 peaktime3 peaktime4 peaktime5 
					pred1 pred2 pred3 pred4 pred5 predtime;
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
		delete afterblackout chmpeak imgpeak integrate maxpct peaks phase2data secpeak peak1-peak5 maxpct
			   phase1calc phase1results phase2results;
	quit; 
%mend recvdwbreast;
*------end of: macro 2: rec.vdw.breast-------------------------;



*------start of: macro 3: rec.vdw--------------------------;
%macro breastrecvdw(inputdata, disease, ptid, month, 
				   sec, chm, rad, hospice, mast, image, obs, cutoffph1=., );
	%if &disease.='breast' %then %do;
		%if &cutoffph1.=. %then %do;
			%let cutoffph1=0.341748894;
		%end;
	%recvdwbreast(inputdata=&inputdata., ptid=&ptid., month=&month., 
				   sec=&sec., chm=&chm., rad=&rad., hospice=&hospice., 
				   mast=&mast., image=&image., obs=&obs., cutoffph1=&cutoffph1.);
	%end;

%mend breastrecvdw;
*--------end of: macro 3: rec.vdw--------------------------;

