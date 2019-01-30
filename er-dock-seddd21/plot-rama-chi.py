#!/usr/bin/python

import sys
from math import *
#import scipy
#import matlot
#import numpy as np
from operator import itemgetter
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
import os

matplotlib.rcParams.update({'font.size': 22})

titles=True
if ((len(sys.argv)>1) and (sys.argv[1]=='notitles')):
        titles=False



def read_data_append(fname,phi,psi,chi1,chi2,phi_ref,psi_ref,chi1_ref,chi2_ref):
	input=open(fname,'r')
	for line in input:
		words=line.split()
		#residue aa_name phi psi chi_1 chi_2 phi psi chi_1 chi_2
		#0       1       2   3   4     5     6   7   8     9
		#first is reference, seocnd is actual
		if ((phi_ref is not None) and (psi_ref is not None)):
			phi_ref.append(float(words[2]))
			psi_ref.append(float(words[3]))
                #it uses "x" to indicate nonexistent chi angles
		if ((chi1_ref is not None) and (chi2_ref is not None) and (words[4]!='x') and (words[5]!='x')):
			chi1_ref.append(float(words[4]))
			chi2_ref.append(float(words[5]))
		if ((phi is not None) and (psi is not None)):
			phi.append(float(words[6]))
			psi.append(float(words[7]))
		#it uses "x" to indicate nonexistent chi angles
		if ((chi1 is not None) and (chi2 is not None) and (words[8]!='x') and (words[9]!='x')):
			chi1.append(float(words[8]))
			chi2.append(float(words[9]))
	input.close()


#def read_data_append_ref(fname,phi,psi,chi1,chi2):
#        input=open(fname,'r')
#        for line in input:
#                words=line.split()
                #residue aa_name phi psi chi_1 chi_2 phi psi chi_1 chi_2
                #0       1       2   3   4     5     6   7   8     9
                #first is starting, seocnd is actual
#                phi.append(float(words[2]))
#                psi.append(float(words[3]))
                #it uses "x" to indicate nonexistent chi angles
#                if ((words[4]!='x') and (words[5]!='x')):
#                        chi1.append(float(words[8]))
#                        chi2.append(float(words[9]))
#        input.close()



#this is the list of target-drug-ref combos to do
l=[('active','gen','1x7r'),('inactive','ral','2qxs')]

states=['hie219','hip219']
script='movemix11-noncmc'
dir=os.environ['HOME']+'/estrogen-receptor2/er-dock-seddd18'
n=60 #number of simulations per state
fnamefmt='{0}/an/rama-{1}-{2}-{3}-{4}-{5}-cutoff3allref-{6}'

ticks=range(-180,181,60)
wt=0.1
path=[(-wt,1),(wt,1),(wt,wt),(1,wt),(1,-wt),(wt,-wt),(wt,-1),(-wt,-1),(-wt,-wt),(-1,-wt),(-1,wt),(-wt,wt)]


for t in l:
	target=t[0]
	drug=t[1]
	ref=t[2]
	print target,drug,ref
	phi=[]
	psi=[]
	chi1=[]
	chi2=[]
	chi1_sc=[]
	chi2_sc=[]
	for state in states:
		for id in range(1,n+1):
			fname=fnamefmt.format(dir,target,state,drug,ref,id,script)
			#print fname
			read_data_append(fname,phi,psi,chi1,chi2,None,None,None,None)
			fname=fnamefmt.format(dir,target,state,drug,ref,id,'sidechainonly')
			#print fname
			read_data_append(fname,None,None,chi1_sc,chi2_sc,None,None,None,None)
			#print phi[-1]
			#print psi[-1]
	phi_ref=[]
        psi_ref=[]
        chi1_ref=[]
        chi2_ref=[]
	fname=fnamefmt.format(dir,target,'hie219',drug,ref,1,script)
	read_data_append(fname,None,None,None,None,phi_ref,psi_ref,chi1_ref,chi2_ref)
	phi_start=[]
	psi_start=[]
	chi1_start=[]
	chi2_start=[]
	if (target=='active'):
		fname=fnamefmt.format(dir,target,'hie219','est','1qku',1,script)
	elif (target=='inactive'):
		fname=fnamefmt.format(dir,target,'hie219','oht','3ert',1,script)
	else:
		print "error"
		sys.exit(-1)
	read_data_append(fname,None,None,None,None,phi_start,psi_start,chi1_start,chi2_start)

	#print sorted(phi)
	#print sorted(psi)
	rama_fig=plt.figure()
	rama_ax=rama_fig.add_subplot('111')
	if (titles):
		rama_ax.set_xlabel(r'$\phi$ (${}^\circ$)')
		rama_ax.set_ylabel(r'$\psi$ (${}^\circ$)')
	rama_ax.set_xlim(left=-180,right=180)
	rama_ax.set_ylim(bottom=-180,top=180)
	rama_ax.set_xticks(ticks)
	rama_ax.set_yticks(ticks)
	rama_ax.grid(b=True,color='k',linestyle='-',linewidth=0.25)
	rama_ax.plot(phi,psi,'g.',label='fully flexible MC',markersize=4)
	rama_ax.plot(phi_ref,psi_ref,'k',label='reference structure',marker=path,linestyle='None',markersize=15,markeredgecolor='None')
	rama_ax.plot(phi_start,psi_start,'r',label='starting structure',marker=path,linestyle='None',markersize=15,markeredgecolor='None')
	rama_ax.legend(fontsize=14,loc='upper right')
	rama_fname='plots2/rama-{0}-{1}-{2}.png'.format(target,drug,ref)
	rama_fig.savefig(rama_fname,output='png',bbox_inches='tight')
	plt.close(rama_fig)
	chi_fig=plt.figure()
	chi_ax=chi_fig.add_subplot('111')
	if (titles):
		chi_ax.set_xlabel(r'$\chi_1$ (${}^\circ$)')
		chi_ax.set_ylabel(r'$\chi_2$ (${}^\circ$)')
	chi_ax.set_xlim(left=-180,right=180)
	chi_ax.set_ylim(bottom=-180,top=180)
        chi_ax.set_xticks(ticks)
        chi_ax.set_yticks(ticks)
        chi_ax.grid(b=True,color='k',linestyle='-',linewidth=0.25)
	chi_ax.plot(chi1,chi2,'g.',label='fully flexible MC',markersize=2)
	chi_ax.plot(chi1_sc,chi2_sc,'b.',label='side chain only',markersize=2)
	chi_ax.plot(chi1_ref,chi2_ref,'k',label='reference structure',marker=path,linestyle='None',markersize=15,markeredgecolor='None')
	chi_ax.plot(chi1_start,chi2_start,'r',label='starting structure',marker=path,linestyle='None',markersize=15,markeredgecolor='None')
	chi_ax.legend(fontsize=14,loc='upper right')
	chi_fname='plots2/chi-{0}-{1}-{2}.png'.format(target,drug,ref)
	chi_fig.savefig(chi_fname,output='png',bbox_inches='tight')
	plt.close(chi_fig)
