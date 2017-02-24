#! /usr/bin/python3.4

import os
import sys
import time
import argparse
import tempfile
import numpy as np
import networkx as nx
from Bio.Align.Applications import MuscleCommandline
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

##########################################################

def blast(seq_file, db_file, proc=4, word_size=10, output_tmpblastfile="", rmout=True, keep_db=False, verbose=2):	
	run=0
	if output_tmpblastfile == "":
		run=1
		output_tmpblastfile=tempfile.NamedTemporaryFile(delete=False).name
		
	output_tmplog=tempfile.NamedTemporaryFile(delete=False).name
	if not os.path.exists(seq_file):
    		return -1
	if not os.path.exists(db_file):
    		return -2
	if rmout == False and verbose > 1:	
		print("     files are not deleted ...")
		print("        Blast output: "+ output_tmpblastfile)
		print("        Blast log: "   + output_tmplog)

	if run == 1 :
		os.system("makeblastdb  -logfile /dev/null -in "+db_file+" -dbtype nucl -parse_seqids")
		
		blastcmd='blastn '
		blastcmd=blastcmd+' '
#		blastcmd=blastcmd+' -num_descriptions 10000000' 
#		blastcmd=blastcmd+' -num_alignments 10000000'
		blastcmd=blastcmd+' -max_target_seqs 10000000'
		blastcmd=blastcmd+' -outfmt 7'
		blastcmd=blastcmd+' -num_threads '   +str(proc)
		blastcmd=blastcmd+' -word_size '   +str(word_size)
		blastcmd=blastcmd+' -query '   +seq_file
		blastcmd=blastcmd+' -db '   +db_file
		blastcmd=blastcmd+' -out '   +output_tmpblastfile
		blastcmd=blastcmd+' -dust yes '
		blastcmd=blastcmd+' 1>2 > '+output_tmplog
		os.system(blastcmd)
		dt=np.dtype([('monomer',np.str,100),('seq',np.str,100),('simil',np.float), ('length',np.int),('x',np.int), ('y',np.int), ('mono_beg',np.int), ('mono_end',np.int), ('seq_beg', np.int), ('seq_end',np.int), ('pvalue', np.float), ('score', np.float)])
		if os.path.getsize(output_tmpblastfile) > 0:
			data = np.genfromtxt(output_tmpblastfile, delimiter='	',
 #    			 names=['monomer','seq','simil', 'length','x', 'y', 'mono_beg', 'mono_end', 'seq_beg', 'seq_end', 'pvalue', 'score']  , 
     			 #dtype=None)
   #  			 dtype=(str,      str,    float, int,      int, int, int,          int,        int,      int,        float,   float))
				 dtype=dt)
			data = np.asarray(data)
		else:
			data=np.array([])
#	print data
	if len(np.shape(data)) == 0:
		data=np.asarray([data], dtype=dt)
	if rmout == True:
		os.remove(output_tmpblastfile)
		os.remove(output_tmplog)
	if keep_db == False:
		os.remove(db_file+".nhr")
		os.remove(db_file+".nin")
		os.remove(db_file+".nsq")
		os.remove(db_file+".nsd")	
		os.remove(db_file+".nog")			
	return(data)
###################################################################################					
def extract_sequences_to_new_file(sequences, db, output_file, mode="w"):
		if mode != "a":
			if os.path.exists(output_file):
				os.remove(output_file)
		for seq in sequences:
			blastcmd='blastdbcmd '
			blastcmd=blastcmd+' -db ' + db
			blastcmd=blastcmd+' -entry ' + seq
			blastcmd=blastcmd+' -dbtype nucl'
			blastcmd=blastcmd+' -outfmt %f'
			blastcmd=blastcmd+' >> '+ output_file
			os.system(blastcmd)
		return(0)

###################################################################################					
def blast_parser(data, min_similarity, min_length):
	subdata = []
	for ii in range(0, len(data)):
		if data[ii][0] != data[ii][1] and data[ii][2] > min_similarity and data[ii][3] > min_length:
			subdata.append([data[ii][0],data[ii][1]])		
	return(subdata)
###################################################################################					
def clusterize(data_pairs):
	clusters = []
	for ii in range(0, len(data_pairs)):
		not_found = True
		for jj in range(0, len(clusters)):
			if data_pairs[ii][0] in clusters[jj] and data_pairs[ii][1] not in clusters[jj]:
				clusters[jj].append(data_pairs[ii][1])
				not_found = False
			elif data_pairs[ii][0] not in clusters[jj] and data_pairs[ii][1] in clusters[jj]:
				clusters[jj].append(data_pairs[ii][0])
				not_found = False
		if not_found == True:
			clusters.append([data_pairs[ii][0], data_pairs[ii][1]])
		
#	for jj in range(0, len(clusters)):
#		clusters[jj]=set(clusters[jj])
	return(clusters)
###################################################################################					
def print_cluster_stats(clusters, verbose=5):
	if verbose > -1:
		print("       Number of clusters: %d" %(len(clusters)))
	if verbose > 0:
		for ii in range(0, len(clusters)):
			print("          Size of cluster %d: %d" % (ii, len(clusters[ii])))	
			if verbose > 3:
				print(clusters[ii])
	return(0)
	
	
###################################################################################						
def clusterize_from_graph(pairs, verbose=3):
		clusters =  []
		G = nx.Graph()
		G.add_edges_from(pairs)
		for cl in nx.connected_components(G):
			clusters.append(list(cl))
		if verbose > 2:
			print("    -> %d clusters" % (len(clusters)))
		return(clusters)
###################################################################################								
def build_all_consensus(clusters, seq_file):
	all_consensus = []
	for ii in range(0, len(clusters)):
		all_consensus.append([build_consensus(clusters[ii], seq_file), build_consensus_name(clusters[ii])])
	return(all_consensus)
###################################################################################								
def build_consensus(cluster, seq_file):
	tmp_file_in = tempfile.NamedTemporaryFile(delete=False).name
	tmp_file_out = tempfile.NamedTemporaryFile(delete=False).name
	
	extract_sequences_to_new_file(cluster, seq_file, tmp_file_in, mode="w")
	
	
	muscle_cline = MuscleCommandline(input=tmp_file_in, out=tmp_file_out,clwstrict=False,quiet=True)

	code = subprocess.call(str(muscle_cline), shell=(sys.platform!="win32"))
	
	alignment = AlignIO.read(tmp_file_out, "fasta")
	summary_align = AlignInfo.SummaryInfo(alignment)
	cons = summary_align.dumb_consensus()

	os.remove(tmp_file_in)
	os.remove(tmp_file_out)
	return(str(cons))
###################################################################################									
def get_len_list_of_lists(mylist):
	l = 0
	for ii in range(len(mylist)):
		l = l + len(mylist[ii])
	return(l)
###################################################################################									
def fusion_list_of_lists(mylist):
	l = []
	for ii in range(len(mylist)):
		for jj in range(len(mylist[ii])):
			l.append(mylist[ii][jj])
	return(l)

	
		
###################################################################################									
def write_unclustered_sequences(seq_file, clusters, output_file, mode="w"):
	clustered_sequences = fusion_list_of_lists(clusters)
	unclustered_sequences = []
	total_sequences = 0
	if mode == "w":
		if os.path.exists(output_file):
			os.remove(output_file)
		
	seq_file_hd = open(seq_file, "r")
	for record in SeqIO.parse(seq_file_hd, "fasta"):
		total_sequences = total_sequences + 1
		if record.id not in clustered_sequences:
			unclustered_sequences.append(record.id)
	extract_sequences_to_new_file(unclustered_sequences, seq_file, output_file, mode = "a")
	return((unclustered_sequences,clustered_sequences))
	
###################################################################################									
def write_all_sequences(seq_file, clusters, all_consensus, output_file, verbose=5):
	(unclust, clust)  = write_unclustered_sequences(args.seq_file, clusters, args.out_file, mode="w")
	write_all_consensus_sequences(all_consensus, args.out_file, mode="a")
	if verbose > 2:
		print("       -> %d unclustered sequences" % (len(unclust)))
		print("       -> %d clustered sequences"   % (len(clust)))
		print("       -> %d consensus sequences"   % (len(all_consensus)))
		print("       -> %d output sequences"      % (len(all_consensus)+len(unclust)))

###################################################################################									
# NOT TESTED
###################################################################################									
def write_clustered_sequences(seq_file, clusters, output_file, mode="w"): 
	clustered_sequences = fusion_list_of_lists(clusters)
	extract_sequences_to_new_file(clustered_sequences, seq_file, output_file, mode = mode)
###################################################################################									
def write_all_consensus_sequences(all_consensus, output_file, mode="w"):
	seq_file_hd = open(output_file, mode)
	for seq,name in all_consensus:
		SeqIO.write(SeqRecord(Seq(seq), id=name), seq_file_hd, "fasta")
	seq_file_hd.close()

###################################################################################								
def build_consensus_name(cluster):
	return(cluster[0] + "_cons_" + str(len(cluster)))
###################################################################################					
parser = argparse.ArgumentParser(description='Build multiple consensus sequences from a fasta file ...', epilog='')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-S', dest="similarity", type=float, default=95.0)
parser.add_argument('-L', dest="length", type=int, default=100)
parser.add_argument('-R', dest="rm_tmp", action='store_false', default=True)
parser.add_argument('-W', dest="word_size", type=int, default=20)
parser.add_argument('-P', dest="processors", type=int, default=4)
parser.add_argument('-o', dest="out_file", required=True)
parser.add_argument('-v', dest="verbose",   type=int, default=0)

args = parser.parse_args()

if not os.path.exists(args.seq_file):
	print("Error: sequence file ("+args.seq_file+") not found !!")
	quit()	

if args.verbose > 0:
	print(" ---> Blasting ...")
data = blast(args.seq_file,args.seq_file, proc=args.processors, word_size=args.word_size, rmout=args.rm_tmp, keep_db=True, verbose=args.verbose)

if args.verbose > 1:
	dim = len(data)
	print("     -> %d hits" % (dim))

if args.verbose > 0:
	print(" ---> Parsing ...")
data_pairs = blast_parser(data, args.similarity, args.length)
if args.verbose > 1:
	dim = len(data_pairs)
	print("     -> %d selected pairs" % (dim))

if args.verbose > 0:
	print(" ---> Clustering ...")
#clusters = clusterize(data_pairs)
clusters  = clusterize_from_graph(data_pairs, args.verbose)
	
if args.verbose > 0:
	print(" ---> Building consensus ...")
all_consensus = build_all_consensus(clusters, args.seq_file)

if args.verbose > 0:
	print(" ---> Writing sequences ...")
write_all_sequences(args.seq_file, clusters, all_consensus, args.out_file)

if args.verbose > 1:
	print(" ---> Printing statistics...")
	print_cluster_stats(clusters, args.verbose)
if args.verbose > 0:
	print(" ---> Bye bye !!!")

