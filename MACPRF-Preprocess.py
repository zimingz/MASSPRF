import numpy as np
import os
import glob
import re
import csv




div_filenames = []
pol_filenames = []
div_lengths = []
sorted_seq = []
pol_replace = []
yeast_data = [[[[[[[]]]]]]]
yeast_data[0].extend(['Pol Gene File', 'Div Gene File', 'Gene length', 'Replacements Pol', 'Replacements Div', 'Synomynous Pol', 'Synomynous Div'])
gene_names = []

#Read in Polymorphism and Divergence Files and get total number of genes
div = glob.glob('./Div*.fas') #creates an array of full filenames 
pol = glob.glob('./Pol*.fas')
size_div = np.size(div)
size_pol = np.size(pol)

#check that two equal for errors:
if size_pol != size_div: print "Error, size do not match for genes"

#Create an Array of Gene Names
for i in range(0,size_pol):
	div_name = os.path.basename(div[i]) #shortens to gene-name.fas
	pol_name = os.path.basename(pol[i])
	div_filenames.append(div_name)
	pol_filenames.append(pol_name)
	np.sort(div_filenames)
	np.sort(pol_filenames)
	
#Create an array of only gene names (i.e. no more .fas or .txt)
for i in range(0,size_pol):
	gene_names.append(os.path.splitext(pol_filenames[i])[0])
	
	
#Find Sizes of Genes and Add to Gene Data Matrix
for i in range(0,size_pol):
	file = open(div_filenames[i], 'r')
	file_data = file.read()
	nucleotide_sequence = re.findall('[ATCG]', file_data)
	length = len(nucleotide_sequence)
	div_lengths.append(length)
	yeast_data.append([pol_filenames[i],div_filenames[i],div_lengths[i],0,0,0,0]) #zeros will be replaced with actual values later
	file.close()
	


#Pre-process Data to get Text Files that List Replacement/Syn Mutations
for i in range(0,size_pol):
	print pol_filenames[i]
	print div_filenames[i]
	os.system('"MAC_PRF_prepv2 -p %s -d %s -o 1 -ci_m 1 -s 1 -t 5 >%s_MACPRF_prep.txt"' %(pol_filenames[i],div_filenames[i],gene_names[i])) #be sure to change -t to est. divergence time for species you're using
	#if no longer fas file, be sure to include species number 

#Get Number of Replacements/Syn
for i in range(0,size_pol):
	with open('%s_MACPRF_prep.txt' %gene_names[i],'r') as file: 
		replace_data = []
		replace_data = ''.join(file.readlines())    
		replace_seq_pol = replace_data[replace_data.index('Polymorphism:'):replace_data.index('Divergence:')]
		replace_seq_div = replace_data[replace_data.index('Divergence:'):replace_data.index('Mission')]
		replacements_pol = replace_seq_pol.count('R')
		replacements_div = replace_seq_div.count('R')
		synomynous_pol = replace_seq_pol.count('S')
		synomynous_div = replace_seq_div.count('S')
		yeast_data[i+1][3] = replacements_pol #Note file info in yeast data starts with row 1 (row zero is column names) 
		yeast_data[i+1][4] = replacements_div
		yeast_data[i+1][5] = synomynous_pol
		yeast_data[i+1][6] = synomynous_div
		file.close()

	

#Create CSV with Gene Info
with open('Yeast_Gene_Info_SynRep.csv', "wb") as csvfile:
	write = csv.writer(csvfile, delimiter = ',')
	for i in range(0,size_div+1):
		write.writerow(yeast_data[i])
	
	
csvfile.close()




	




 


	

	

