#!/usr/local/bin/python
from __future__ import division
import syssys.path.insert(0, "/Applications/biopython-1.63")
import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid
import csvimport sys
import re
from re import *
from Bio.Seq import Seq
from StringIO import StringIO
import itertools
from itertools import combinations
from Bio.Alphabet import IUPAC
import os
import subprocess
import time
import matplotlib.pyplot as plt
import svgutils.transform as sg
import math
import glob
from collections import defaultdict
import operator
# new blast-free version of SLIM, uses ublast and protein database
# SLIM_ubp3 uses full path to data base
# SLIM_ubp4 adds steps to make hit file unique and then selects only best hits

start1 = time.time()
if len(sys.argv)!=4:
	print("Usage: python SLIM_ubp.py contigs.fasta data_base_list filter500_yes_no")
	sys.exit()	
print "Running SLIM_ubp.py"

input_file = sys.argv[1]
summary_outprefix = os.path.splitext(input_file)[0]
outprefix = os.path.splitext(input_file)[0]
data_base_list_file= sys.argv[2]
data_base_list = os.path.splitext(data_base_list_file)[0]
turn_on_size_filter=sys.argv[3] 
print "counting and making you wait"
records = list(SeqIO.parse(input_file, "fasta"))number_contigs = len(records)

unique_hit_fasta_list=[] #to collect files for graphing contig lengths
virus_family_names=[]
virus_genome_lengths=[]

#Set of functions to generate rev_com (without using biopython/fasta)
def reverse(s): 
	letters = list(s) 
	letters.reverse() 
	return ''.join(letters) 
def complement(s): 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y', 'K':'M','M':'K', 'S':'S', 'W':'W', 'N':'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'y':'r', 'r':'y', 'k':'m','m':'k', 's':'s', 'w':'w', 'n':'n'} 
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters) 
def rev_com(s): 
	s = reverse(s) 
	s = complement(s) 
	return s

print "printing details to summary file"
summary_file = open(summary_outprefix+'_mapped_to_'+data_base_list+'_summary.txt','a')
print >> summary_file, 'Contig_file: %s\nNumber_contigs: %i\nVirus_data_base: %s' %(input_file, number_contigs,data_base_list_file)
summary_file.close()

#add size cutoff at 500nt and above if number of contigs >3000
if turn_on_size_filter == "yes":
	if number_contigs >3000:
		call_string_rmv_shorties = "python Remove_shorties_sz.py "+input_file+" 499" 
		subprocess.call(call_string_rmv_shorties, shell=True)
		input_file = outprefix+"_gr_499.fas"
		outprefix = os.path.splitext(input_file)[0]
		records_size_filtered = list(SeqIO.parse(input_file, "fasta"))		number_contigs_gr499 = len(records_size_filtered)
		summary_file = open(summary_outprefix+'_mapped_to_'+data_base_list+'_summary.txt','a')
		print >> summary_file, 'Clustered_Contig_file: %s\nNumber_contigs: %i' %(input_file, number_contigs_gr499)
		summary_file.close()

#generate virus_database_list
virus_database_list =[]
with open(data_base_list_file, "rU") as f:
	reader = csv.reader(f, quoting=csv.QUOTE_NONE)
#	reader.next()#skip header
	for row in reader:			
		virus_database_list.append(row[0])
		virus_family_names.append(row[1])
		virus_genome_lengths.append(row[2])		
f.close()
print virus_database_list
print virus_family_names
print virus_genome_lengths

summary_file = open(summary_outprefix+'_mapped_to_'+data_base_list+'_summary.txt','a')
print >> summary_file, 'virus_database_list: %s' %virus_database_list
summary_file.close()

#number_virus_families_searched = len(virus_database_list)

for vdb in virus_database_list:
	vdb_path_parts = vdb.split('/')
	vdb_parts = vdb_path_parts[-1].split('.')
	vdb_outprefix = vdb_parts[0]
	
	print "Searching in "+vdb_path_parts[-1]
		
	ublast_call_string ="/Applications/usearch -ublast "+input_file+" -db "+vdb+" -evalue 1e-9 -accel 0.5 -alnout "+outprefix+"_"+vdb_outprefix+"_hits_aln.txt -matched "+outprefix+"_"+vdb_outprefix+"_hits_matched.txt -userout "+outprefix+"_"+vdb_outprefix+"_simple.txt -userfields query+qlo+qhi+id+qframe+target -maxhits 1"

	subprocess.call(ublast_call_string, shell=True)
	
	ublast_output_table = outprefix+"_"+vdb_outprefix+"_simple.txt"
	
	summary_file = open(summary_outprefix+'_mapped_to_'+data_base_list+'_summary.txt','a')
	print >> summary_file, 'usearch_output_table: %s' %ublast_output_table
	summary_file.close()
	#parse ublast_output_table
	contig_hit_list =[]
	ublast_reader = csv.reader(open(ublast_output_table, "rU"), delimiter='\t')
	for row in ublast_reader:
		contig_hit_list.append(row)

	print "sorting contigs with " +vdb_outprefix+ " hits"
	hit_fasta = open(outprefix+"_"+vdb_outprefix+"_hits.fas", "a")	
	hit_fasta.close() #need to make file even if it stays empty
	hit_fasta = open(outprefix+"_"+vdb_outprefix+"_hits.fas", "a")	
	for record in SeqIO.parse(open(input_file, "rU"), "fasta"):
		for element in contig_hit_list:
			if record.id == element[0]:	
				nopipe_name = element[0].replace("|", "_")		
				if int(element[4]) > 0:
					new_name = nopipe_name+"_"+element[5]
					print >> hit_fasta, '>'+new_name
					print >> hit_fasta, record.seq
				elif int(element[4]) <0:
					new_name_rc = nopipe_name+"_rc_"+element[5]
					rev_com_seq = rev_com(record.seq)
					print >> hit_fasta, '>'+new_name_rc
					print >> hit_fasta, rev_com_seq
	hit_fasta.close() 


	print "preparing unique set of contigs with " +vdb_outprefix+ " hits"
	harvested_list =[]
	f2 = open(outprefix+"_"+vdb_outprefix+"_unique_hits.fas", "a")
	f2.close()#need to make file even if it stays empty
	for record in SeqIO.parse(open(outprefix+"_"+vdb_outprefix+"_hits.fas", "rU"), "fasta"):
		record_name = record.id
		name_parts = record_name.split('_')
		if len(name_parts)>=4:
			contig_name = name_parts[0]+"_"+name_parts[1]+"_"+name_parts[2]+"_"+name_parts[3]
		else:
			contig_name = record_name
		if contig_name not in harvested_list: 
			harvested_list.append(contig_name)
	 		f2 = open(outprefix+"_"+vdb_outprefix+"_unique_hits.fas", "a")
			print >> f2, '>'+'%s' % record_name
			print >> f2, record.seq
			f2.close()
	unique_hit_fasta_list.append(outprefix+"_"+vdb_outprefix+"_unique_hits.fas")

#collect contigs lengths
contig_length_list_collector=[]
for i in range(len(unique_hit_fasta_list)):
	contig_length_list=[]
	contig_file = unique_hit_fasta_list[i] 
	outprefix2 = os.path.splitext(contig_file)[0]
	virus_screened= virus_family_names[i]
	genome_length= virus_genome_lengths[i]
	if os.path.getsize(unique_hit_fasta_list[i]) == 0:#check if file has contents
		print virus_screened +": empty file no contigs"
		contig_length_list.append(0)#need to make list for empty files		
	else:
		for record in SeqIO.parse(open(contig_file, "rU"), "fasta"):
			contig_length = len(record.seq)
			contig_length_list.append(contig_length)
	contig_length_list_collector.append((contig_length_list))
	
print contig_length_list_collector

#add graphing 
fig = plt.figure()
index=1
page_number=0
total_counter=0
for i in range(len(contig_length_list_collector)):
	contig_file = unique_hit_fasta_list[i] 
	outprefix2 = os.path.splitext(contig_file)[0]
	virus_screened= virus_family_names[i]
	genome_length= virus_genome_lengths[i]				
	total_number_contigs = len(contig_length_list_collector[i])
	largest_contig = max(contig_length_list_collector[i])
	print "total_number_contigs: "+str(total_number_contigs)
	print "largest contig: " +str(largest_contig)
	ax2 = fig.add_subplot(6,1,index)
#	if largest_contig == 0:
#		number_of_bins= int(int(genome_length)/30)
#	else:
#		number_of_bins=int(int(largest_contig)/30)
	number_of_bins = 3000
	ax2.hist(contig_length_list_collector[i], bins = number_of_bins, color='#36486b', edgecolor = '#36486b') #darkblue reduce bin number to reduce figure size		

#	ax2.hist(contig_length_list, bins = int(largest_contig), color='#36486b', edgecolor = '#36486b') #darkblue		
	ax2.set_yticklabels(["","","1","","2","","3","",">=4"])
	ax2.tick_params(axis='both', which='major', labelsize=8)		
	ax2.patch.set_facecolor('#d1e0e0')#a green grey
	ax2.patch.set_alpha(0.6)
	ax2.set_xlabel(virus_screened+"_"+outprefix+"_contig_length", fontsize=8)
	ax2.set_ylabel("number of contigs", fontsize=6)
	plt.ylim(0, 4)
	if largest_contig == 0:
		upper_limit= int(genome_length)+100
	else:
		upper_limit = largest_contig+100
	plt.xlim(0, upper_limit)				
	almost_full_genome = int(0.9*(int(genome_length)))#add line marking 90% full genome size
	plt.axvline(almost_full_genome, color='red', linestyle='dotted', linewidth=2)
	
	index = index+1
	total_counter = total_counter+1
	if index==7 or total_counter == len(contig_length_list_collector):
		page_number = page_number+1
		plt.tight_layout()
		plt.savefig(outprefix+"_histogram_contigs_pg"+str(page_number)+".pdf")
		fig = plt.figure()		
		index=1
	
		
print "That's All Folks!"
		
end1 = time.time()
elapsed1 = round(end1 - start1, 2)elapsed1min = round(elapsed1/60, 2)print "Elapsed time (minutes): ",elapsed1min

summary_file = open(summary_outprefix+'_mapped_to_'+data_base_list+'_summary.txt','a')
print >> summary_file, "Elapsed time (minutes): ",elapsed1min
summary_file.close()	
	
print 'That\'s All Folks!'