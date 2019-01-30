#!/usr/bin/python

import sys
from math import *
#import scipy
#import matlot
import numpy as np
from operator import itemgetter
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
import os

def get_y(l):
        n=len(l)
        return [i/float(n) for i in range(1,n+1)]

matplotlib.rcParams.update({'font.size': 22})

titles=True
if ((len(sys.argv)>1) and (sys.argv[1]=='notitles')):
	titles=False

#ctype=sys.argv[1]

extpath='/net/roos/home4/zuckerman/spiriti/estrogen-receptor2/er-dock-seddd18'
scripts=['cutoff3allref-movemix11a','cutoff3allref-movemix11-noncmc','cutoff3allref-proteinfixed','cutoff3allref-sidechainonly']
dir = {}
dir['cutoff3allref-movemix11a']='.'
dir['cutoff3allref-movemix11-noncmc']=extpath
dir['cutoff3allref-proteinfixed']=extpath
dir['cutoff3allref-sidechainonly']=extpath
maxrank=10 #for figs 8(a) and (c)
ranks=range(1,maxrank+1)

druglistfname="matched-working-drugs"
druglistfile=open(druglistfname,'r')
drugs = []
targets = {}
#refs = {}
all_drugrefs = [] 
drugref_by_target = {}
for drugline in druglistfile:
	words=drugline.split()
	target=words[0]
	drug=words[1]
	ref=words[2]
	if (drug not in drugs):
		drugs.append(drug)
	targets[drug]=target
	if (target not in drugref_by_target):
		drugref_by_target[target]=[]
	t=(drug,ref)
	drugref_by_target[target].append(t)
	if (t not in all_drugrefs):
		all_drugrefs.append(t)
	#refs[drug]=ref
	
druglistfile.close()

#print drugs
#print targets
#drugref_by_target = {target: [] for target in targets}
#print drugref_by_target
#quit()


ctypes = ["3a", "0.25p", "0.10p"]
enames = {3: 'intxn-energy', 4:'total-energy'}
bestrmsds_lg_cluster_by_script={}
bestrmsds_le_cluster_by_script={}
bestrmsds_no_cluster={}
avg_rmsd_rank_by_script={}
target_labels={'active':'agonists', 'inactive':'antagonists'}


for ctype in ctypes:
	sims_by_drugref = {}
	rmsds = {}
	#print drugrefs
	for script in scripts:
		for t1 in all_drugrefs:
			(drug, ref) = t1
			target = targets[drug]
			try:
				summaryfname='{0}/an/summary-{1}-{2}-{3}-{4}'.format(dir[script],target,drug,ref,script)
				summaryfile=open(summaryfname,'r')
				clusterfname='{0}/an/clusters-drugrmsd-{1}-{2}-{3}-{4}'.format(dir[script],ctype,target,drug,script)
				clusterfile=open(clusterfname,'r')
			except:
				print "data missing for ",target,drug,ref
				continue
			#os.system("wc -l {0}".format(summaryfname))
			#os.system("wc -l {0}".format(clusterfname))
			for (sline, cline) in zip(summaryfile,clusterfile):
				sumwords=sline.split()
				cluswords=cline.split()
				try:
					#these field numbers are zero-based, corresponding to fields 2 (state), 5 (id), 7 (RMSD), 18 
					state=sumwords[1]
					id=int(sumwords[4])
					rmsd=float(sumwords[6])
					intxn_energy=float(sumwords[17])
					total_energy=float(sumwords[13])
					cluster=int(cluswords[1])
				except:
					print summaryfname, sline
					continue

				t2 = (drug,ref,script)
				if (t2 not in sims_by_drugref):
					sims_by_drugref[t2]=[]
                        	#   0     1  2    3            4            5
				t3=(state,id,rmsd,intxn_energy,total_energy,cluster)
				sims_by_drugref[t2].append(t3)
			summaryfile.close()
			clusterfile.close()
			if (len(sims_by_drugref[t2])<120):
				print "warning: failed to read data from ",summaryfname

#print sims_by_drugref

#script="cutoff3allref-movemix11a"
	for efld in [3,4]:
		for script in scripts:
			for target in drugref_by_target:
				drugrefs=drugref_by_target[target]
				bestrmsds_min_intxn_energy = []
				bestrmsds_min_intxn_energy_largest_cluster = []
				bestrmsds_min_intxn_energy_lowest_energy_cluster = []
				avg_rmsd_by_rank = [0]*maxrank
				for drugref in drugref_by_target[target]:
					t2=drugref+(script,)
					sims=sims_by_drugref[t2]
					#if (drugref[0]=='llc'):
						#print "###",script,target,len(sims),[(t[2],t[3]) for t in sims]
						#print "###",script,target,[t[3] for t in sims]
					#sort by intxn energy, pick the lowest, obtain the RMSD
					sorted_sims=sorted(sims,key=itemgetter(efld))
					#print sorted_sims
					nclus=max(sims,key=itemgetter(5))[5]
					#print drugref,nclus
					#continue
					bestrmsd_min_intxn_energy = sorted_sims[0][2]
					#print drugref,bestrmsd_min_intxn_energy
					bestrmsds_min_intxn_energy.append(bestrmsd_min_intxn_energy)
					#minimum intxn energy in largest cluster (cluster #1)
					sims_largest_cluster = [t for t in sorted_sims if (t[5]==1)]
					bestrmsd_min_intxn_energy_largest_cluster = sims_largest_cluster[0][2]
					bestrmsds_min_intxn_energy_largest_cluster.append(bestrmsd_min_intxn_energy_largest_cluster)
					#average energy by cluster -- there might be a clerver way to do this in python
					avg_energy_by_cluster = [0.0]*nclus
					count_by_cluster = [0]*nclus
					for t in sims:
						iclus=t[5]
						avg_energy_by_cluster[iclus-1] = avg_energy_by_cluster[iclus-1] + t[efld]
						count_by_cluster[iclus-1] = count_by_cluster[iclus-1]+1
					#print drugref,count_by_cluster
					for iclus in range(0,nclus):
						if (count_by_cluster[iclus]>0):
							avg_energy_by_cluster[iclus]=avg_energy_by_cluster[iclus]/count_by_cluster[iclus]
						else:
							avg_energy_by_cluster[iclus]=1.0e20
					#print drugref, avg_energy_by_cluster
					avg_energy_by_cluster=zip(range(1,nclus+1),avg_energy_by_cluster)
					lowest_energy_cluster=min(avg_energy_by_cluster,key=itemgetter(1))[0]
					#print drugref,lowest_energy_cluster
					sims_lowest_energy_cluster = [t for t in sorted_sims if (t[5]==lowest_energy_cluster)]
					bestrmsd_min_intxn_energy_lowest_energy_cluster = sims_lowest_energy_cluster[0][2]
					bestrmsds_min_intxn_energy_lowest_energy_cluster.append(bestrmsd_min_intxn_energy_lowest_energy_cluster)
					#print drugref, bestrmsd_min_intxn_energy_lowest_energy_cluster
					#if ((drugref[0]=='llc')): # or (drugref==('oht','3ert'))):
						#print "***", script, target, [t[2] for t in sorted_sims]
						#print "***", script, target, [t[3] for t in sorted_sims]
	
					for rank in ranks:
						#select  the N best poses by interaction energy 
						#and select the lowest RMSD from among them.
						#print sorted_sims[0:rank]
						min_rmsd_rank=min(sorted_sims[0:rank],key=itemgetter(2))[2]
						if ((script=="cutoff3allref-movemix11-noncmc") and (target=="inactive")):
							print "***", script, target, drugref, rank, min_rmsd_rank
						#average by rank
						avg_rmsd_by_rank[rank-1] = avg_rmsd_by_rank[rank-1] + min_rmsd_rank
				
					avg_rmsd_by_rank = np.array(avg_rmsd_by_rank)/float(len(drugref_by_target[target]))
					print script,target,avg_rmsd_by_rank
				avg_rmsd_rank_by_script[(target,script)]=avg_rmsd_by_rank

				#to produce cdf plots, must sort each list	
				bestrmsds_min_intxn_energy = sorted(bestrmsds_min_intxn_energy)
				bestrmsds_min_intxn_energy_largest_cluster = sorted(bestrmsds_min_intxn_energy_largest_cluster)
				bestrmsds_min_intxn_energy_lowest_energy_cluster = sorted(bestrmsds_min_intxn_energy_lowest_energy_cluster)

				bestrmsds_lg_cluster_by_script[(ctype,target,script)]=bestrmsds_min_intxn_energy_largest_cluster
				bestrmsds_le_cluster_by_script[(ctype,target,script)]=bestrmsds_min_intxn_energy_lowest_energy_cluster
				bestrmsds_no_cluster[(target,script)]=bestrmsds_min_intxn_energy
				avg_rmsd_rank_by_script[(target,script)]=avg_rmsd_by_rank

				#print target,script,bestrmsds_min_intxn_energy #,get_y(bestrmsds_min_intxn_energy)
				print "making plots for",target,ctype,script
				fig7=plt.figure()
	        		ax7=fig7.add_subplot(111)
				if (titles):
					ax7.set_xlabel("RMSD (A)")
					ax7.set_ylabel("cumulative probability")
	                		ax7.text(0.05,0.95,target_labels[target],fontsize=24,transform=ax7.transAxes,
        	                		horizontalalignment='left',verticalalignment='top')

				ax7.plot(bestrmsds_min_intxn_energy,get_y(bestrmsds_min_intxn_energy),
					'r-',label='min interaction energy overall',linewidth=2)
				ax7.plot(bestrmsds_min_intxn_energy_largest_cluster,get_y(bestrmsds_min_intxn_energy_largest_cluster),
					'g-',label='using largest cluster',linewidth=2)
				ax7.plot(bestrmsds_min_intxn_energy_lowest_energy_cluster,get_y(bestrmsds_min_intxn_energy_lowest_energy_cluster),
					'b-',label='using lowest energy cluster',linewidth=2)
				ax7.legend(fontsize=18,loc='lower right')
				#ax7.legend_.remove()

				plotfname='plots/scoring-methods-{0}-{1}-{2}-{3}.png'.format(enames[efld],target,ctype,script)
				fig7.savefig(plotfname,output='png',bbox_inches='tight')
				plt.close(fig7)

			
quit()
#print bestrmsds_no_cluster
#print bestrmsds_le_cluster_by_script
#print avg_rmsd_rank_by_script
#figure 8(b) and (d): CDF of best rmsds by script
script_labels={"cutoff3allref-movemix11a": "mixed NCMC/MC", "cutoff3allref-movemix11-noncmc": "fully flexible MC",
	"cutoff3allref-proteinfixed": "fixed protein",	"cutoff3allref-sidechainonly": "side chain only"}
script_colors={"cutoff3allref-movemix11a": "m-", "cutoff3allref-movemix11-noncmc": "g-",
	"cutoff3allref-proteinfixed": "r-", "cutoff3allref-sidechainonly": "b-"}

ctype_labels={"3a":"Constant 3 A cluster radius", "0.25p": "25th percentile", "0.10p": "10th percentile"}
ctype_colors={"3a":"r-","0.25p":"g-","0.10p":"b-"}


for target in drugref_by_target:
	for ctype in ctypes:
		fig8bd=plt.figure()
		ax8bd=fig8bd.add_subplot(111)
		ax8bd.set_xlabel("RMSD (A)")
		ax8bd.set_ylabel("cumulative probability")
		ax8bd.text(0.05,0.95,target_labels[target],fontsize=24,transform=ax8bd.transAxes,
			horizontalalignment='left',verticalalignment='top')
		for script in scripts:
			t = (ctype,target, script)
			ax8bd.plot(bestrmsds_le_cluster_by_script[t],get_y(bestrmsds_le_cluster_by_script[t]),
				script_colors[script],label=script_labels[script],linewidth=2)
		ax8bd.legend(fontsize=14,loc='lower right')
		plotfname='plots/cum-dist-best-rmsd-{0}-{1}.png'.format(target,ctype)
		fig8bd.savefig(plotfname,output='png',bbox_inches='tight')
		plt.close(fig8bd)

        fig8bd=plt.figure()
        ax8bd=fig8bd.add_subplot(111)
        ax8bd.set_xlabel("RMSD (A)")
        ax8bd.set_ylabel("cumulative probability")
        ax8bd.text(0.05,0.95,target_labels[target],fontsize=24,transform=ax8bd.transAxes,
                horizontalalignment='left',verticalalignment='top')
        for script in scripts:
                t = (target, script)
                ax8bd.plot(bestrmsds_no_cluster[t],get_y(bestrmsds_no_cluster[t]),
                	script_colors[script],label=script_labels[script],linewidth=2)
        ax8bd.legend(fontsize=14,loc='lower right')
        plotfname='plots/cum-dist-best-rmsd-{0}-nocluster.png'.format(target)
        fig8bd.savefig(plotfname,output='png',bbox_inches='tight')
        plt.close(fig8bd)


	fig8ac=plt.figure()
	ax8ac=fig8ac.add_subplot(111)
	ax8ac.set_xlabel("rank")
	ax8ac.set_ylabel("best RMSD averaged over ligands")
 	ax8ac.text(0.05,0.95,target_labels[target],fontsize=24,transform=ax8ac.transAxes,
		horizontalalignment='left',verticalalignment='top')
	#ax8ac.set_ylim(ymin=0.0,auto=True)

	for script in scripts:
		t=(target,script)
		ax8ac.plot(ranks,avg_rmsd_rank_by_script[t],
			script_colors[script],label=script_labels[script],linewidth=2)
	ax8ac.legend(fontsize=14,loc='upper right')
	ax8ac.set_ylim(bottom=0.0)
	plotfname='plots/avg-best-cum-rmsd-{0}-intxn-energy.png'.format(target)
	fig8ac.savefig(plotfname,output='png',bbox_inches='tight')
	plt.close(fig8ac)


        comp_fig=plt.figure()
        comp_ax=comp_fig.add_subplot(111)
        comp_ax.set_xlabel("RMSD (A)")
        comp_ax.set_ylabel("cumulative probability")
        for ctype in ctypes:
		t = (ctype,target, "cutoff3allref-movemix11a")
		comp_ax.plot(bestrmsds_le_cluster_by_script[t],get_y(bestrmsds_le_cluster_by_script[t]),
			ctype_colors[ctype],label=ctype_labels[ctype],linewidth=2)
	t=(target,"cutoff3allref-movemix11a")
        comp_ax.plot(bestrmsds_no_cluster[t],get_y(bestrmsds_no_cluster[t]),
        	'k-',label="no clustering",linewidth=2)

        comp_ax.legend(fontsize=14,loc='lower right')
        compplotfname='plots/cluster-methods-{0}.png'.format(target)
        comp_fig.savefig(compplotfname,output='png',bbox_inches='tight')
        plt.close(comp_fig)


#print  sims_by_drugref

