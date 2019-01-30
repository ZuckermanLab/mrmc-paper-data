function logsum(a,b) # ln(e^a+e^b)
{
        #print a,b;
        max=a; absdiff=a-b;
        if (b>a) {max=b; absdiff=b-a};
        #print max,absdiff;
        #get rid of annoying "out of range" errors
        if (absdiff>50) return max; else return max+log(1+exp(-absdiff));
	f=cycle/block;
}


BEGIN {
	kt=0.596;
	#if (type~/trans/) angle=0; else angle=1; 
	#if (angle) {nbin=18; dbin=10} else {nbin=50; dbin=0.01};
	for (ibin=0; ibin<nbin; ibin++) {count[ibin]=0; sumprob_ncmc[ibin]=-1000; sumprob_mc[ibin]=-1000};
}

#clear data from one simulation
/END OF SIMULATION/ {delete p_ncmc}

#actually stores logarithms; we are working in the log domain
(NF==2) {
	if ($2<0) p=0.0; else p=-$2/kt;
	p_ncmc[$1]=p;
	#print $1,p;
}

(NF>2) && (tolower($2)==tolower(type)) {
        #identify last step of NCMC cycle involving this step
        step=$1;
	size=$3;
        if ($4<0) p_intramove=0.0; else p_intramove=-$4/kt; #actually log(p)
        x=step/block;
        n=int(x);
	rstep=step%block;
	if (rstep>cycle) {
		#MC-only phase of NCMC cycle
		sumprob_mc[bin]=logsum(sumprob_mc[bin],p_intramove);
		#if ((size>=60) && (size<=105)) print size,rstep,p_intramove;
	} else {
		#NCMC in use
 		ncmc_step=n*block+cycle;
        	#if (frac>0) n++;
       		#ncmc_step=n*cycle;
        	#if ((step%1000)==0) 
		#print step,ncmc_step,p_ncmc[ncmc_step];
        	bin=int(size/dbin);
        	count[bin]++;
        	sumprob_ncmc[bin]=logsum(sumprob_ncmc[bin],p_intramove+p_ncmc[ncmc_step]);
	}
}


END {
	for (ibin=0; ibin<nbin; ibin++) {
		center=(ibin+0.5)*dbin;
		if (count[ibin]>0) {
			avg_ncmc=exp(sumprob_ncmc[ibin]-log(count[ibin])); 
			avg_mc=exp(sumprob_mc[ibin]-log(count[ibin]));
		} else {
			avg_ncmc=0.0;
			avg_mc=0.0;
		}
		print center,count[ibin],avg_ncmc,avg_mc;#exp(sumprob[ibin]);
	}
}
