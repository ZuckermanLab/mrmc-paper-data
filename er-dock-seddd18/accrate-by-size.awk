BEGIN {
	#if (type~/trans/) angle=0; else angle=1; 
	#if (angle) {nbin=18; dbin=10} else {nbin=50; dbin=0.01};
	for (ibin=0; ibin<nbin; ibin++) {count[ibin]=0; sumprob[ibin]=0.0};
	total=0;
	lines=0;
}

(tolower($2)==tolower(type)) {
	size=$3;
	p=$5;
	bin=int(size/dbin);
	count[bin]++;
	sumprob[bin]+=p;
	total++;
}

{lines++} #count all of the moves, regardless of type
#"sumprob" = count*avgprob = expected number of moves of this type and size accepted per simulation.
END {
	cdf=0;
	for (ibin=0; ibin<nbin; ibin++) {
		center=(ibin+0.5)*dbin;
		if (count[ibin]>0) avg=sumprob[ibin]/count[ibin]; else avg=0.0;
		cdf+=count[ibin]/total;
		print center,count[ibin],avg,cdf,sumprob[ibin]/lines;
	}
}
