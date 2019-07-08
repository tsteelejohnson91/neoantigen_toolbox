#!/bin/python
import time
print('Start')
start_time = time.time()
import os
import numpy as np
import sys 
import subprocess as sp
import multiprocessing as mp

def uuid_map(fname, cancer):
	with open(fname) as f:
		lines = f.read().splitlines()
	lines = [line.rstrip('\n').split('\t') for line in lines]
	patients = list()
	for line1 in lines:
		if line1[3] == 'Aligned Reads' and line1[7] == 'Solid Tissue Normal':
			for line2 in lines:
				if line2[3] == 'Aligned Reads' and line2[7] == 'Primary Tumor' and line1[5] == line2[5]:
					patients.append([line1[5], line2[0], line2[1], line1[0], line1[1]])
	fout = open(cancer+'_uuids.txt','w')
	for patient in patients:
		fout.write(patient[0]+'\t'+patient[1]+'\t'+patient[2]+'\t'+patient[3]+'\t'+patient[4]+'\n')
	return('Done')

def bam2seq(args):
	print('Converting bams to sequences')
	fname = args[0]
	oname = args[1]
	print('module load samtools/1.2; samtools view '+fname+' > '+fname[0:fname.find('.')]+'.sam')
	os.system('module load samtools/1.2; samtools view '+fname+' > '+fname[0:fname.find('.')]+'.sam')
	print('awk \'{print $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
	os.system('awk \'{print $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
	print('rm '+fname[0:fname.find('.')]+'.sam')
	os.system('rm '+fname[0:fname.find('.')]+'.sam')
	print('Done: '+str(round(time.time() - start_time,4))+'s')
	os.system('module unload samtools/1.2')
	return('Done')

def fastq2seq(args):
	print('Converting fastqs to sequences')
	fname1 = args[0]
	fname2 = args[1]
	oname = args[2]
	os.system('awk "(NR%4==2)" '+fname1+' > tmp1.txt')
	os.system('awk "(NR%4==2)" '+fname2+' > tmp2.txt')
	os.system('cat tmp1.txt tmp2.txt > '+oname)
	os.system('rm tmp1.txt tmp2.txt')
	return('Done')

def bam2fastq(args):
	# convert bam file to fastqs (MIND-NUMBINGLY SLOW)
	fpath = args[0]
	fname_prefix = args[1]
	line1 = 'module load gcc;'
	line2 = 'module load bedtools;'
	line3 = 'module load samtools/1.2;'
	line4 = 'samtools sort -n '+fpath+'/'+fname_prefix+'.bam '+fpath+'/aln.qsort;'
	line5 = 'bedtools bamtofastq -i '+fpath+'/aln.qsort.bam -fq '+fpath+'/'+fname_prefix+'_rd1.fq -fq2 '+fpath+'/'+fname_prefix+'_rd2.fq;'
	line6 = 'module unload gcc;'
	line7 = 'module unload bedtools;'
	line8 = 'module unload samtools/1.2;'
	os.system(line1+line2+line3+line4+line5+line6+line7+line8)

def fastq2hla(args):
	# Generating HLA predictions
	fpath = args[0]
	fname_prefix = args[1]
	tool_path = args[2]
	pat_num = args[3]
	line1 = 'module unload PrgEnv-cray;'
	line2 = 'module load PrgEnv-intel;'
	line3 = 'module load r;'
	line4 = 'PATH=$PATH:'+tool_path+'/bowtie-0.12.7;'
	line5 = 'PATH=$PATH:'+tool_path+'/seq2HLA;'
	line6 = 'chmod u+x '+tool_path+'/seq2HLA/seq2HLA.py;'
	line7 = 'chmod u+x '+tool_path+'/seq2HLA/command.R;'
	line8 = 'PYTHONPATH=$PYTHONPATH:'+tool_path+'/biopython-1.58;'
	line9 = 'cd '+tool_path+'/seq2HLA; outputPrefix="hlaPred'+str(pat_num)+'";'
	line10 = 'seq2HLA.py -1 '+fpath+'/'+fname_prefix+'_rd1.fq -2 '+fpath+'/'+fname_prefix+'_rd2.fq -r $outputPrefix -l 50;'
	line11 = 'module unload r;'
	line12 = 'module unload PrgEnv-intel;'
	line13 = 'module load PrgEnv-cray;'
	os.system(line1+line2+line3+line4+line5+line6+line7+line8+line9+line10+line11+line12+line13)
	# Extracting HLA predictions from output file
	HLAalleles = []
	fin = open(tool_path+'/seq2HLA/hlaPred'+str(pat_num)+'-ClassI.HLAgenotype','r')
	line = fin.readline()
	line = fin.readline()
	while line != '':
		line = line.rstrip('\n')
		if 'hoz' in line:
			line = line.split('\t')
			tmp = line[1]
			tmp = tmp.split('*')
			HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
		else:
			line = line.split('\t')
			tmp = line[1]
			tmp = tmp.split('*')
			HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
			tmp = line[3]
			tmp = tmp.split('*')
			HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
		line = fin.readline()
	return HLAalleles

def split_file(args):
	fname = args[0]
	nfiles = int(args[1])
	# Begin by counting the lines of input sequence file
	print('Counting lines')
	if os.path.isfile(fname+'.txt'):
		os.system('wc -l '+fname+'.txt > '+fname+'_tmp.txt')
	elif os.path.isfile(fname):
		os.system('wc -l '+fname+' > '+fname+'_tmp.txt')
	else:
		print('Error: function: split_file(): file does not exist')
	with open(fname+'_tmp.txt') as f:
		N = f.readline()
	os.system('rm '+fname+'_tmp.txt')
	N = N.split(" ")
	N = int(N[0])
	print('Done: '+str(round(time.time() - start_time,4))+'s \tLines: '+str(N))
	# Split input sequence file into equal sized subfiles for multiprocessing
	print('Splitting into multiple files')
	if os.path.isfile(fname+'.txt'):
		os.system('split -l '+str(N/nfiles+1)+' '+fname+'.txt '+fname)
	else:
		os.system('split -l '+str(N/nfiles+1)+' '+fname+' '+fname)
	print('Done: '+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
	return('Done')

def merge_files(args):
	fpath = args[0]
	fname_prefix = args[1]
	nfiles = args[2]
	exten = args[3]
	print('Merging ' + fname_prefix + ' files into one file')
	alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
	i = 0;
	j = 0;
	fnames = [''] * nfiles
	#print(fnames)
	while 25*i+j < nfiles-1:
		#print(str(i)+' '+str(j))
		#print(fnames[25*i+j])
		fnames[25*i+j] = fpath+'/'+fname_prefix+alphabet[i]+alphabet[j]+exten
		j = j + 1
		if j == 26:
			j = 0
			i = i+1
	#print(fnames)
	command = 'cat'
	for i in range(nfiles):
		command = command + ' ' + fnames[i]
	command = command + ' > ' + fpath + '/' + fname_prefix + exten
	print(command)
	os.system(command)
	print('Done: '+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
	return('Done')

def translate(seq): 
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 

def translate_file(fname_prefix_sub):
	infile = open(fname_prefix_sub, 'r')
	outfile = open(fname_prefix_sub+'_amino', 'w')
	i=0;
	line = infile.readline()
	line = line.rstrip('\n')
	while line != '':
		if 'N' not in line:
			tmp = line[0:(len(line)-(len(line)%3))]
			aa1 = translate(tmp)
			tmp = line[1:(len(line)-((len(line)-1)%3))]
			aa2 = translate(tmp)
			tmp = line[2:(len(line)-((len(line)-2)%3))]
			aa3 = translate(tmp)
			outfile.write(aa1+'\n')
			outfile.write(aa2+'\n')
			outfile.write(aa3+'\n')
		line = infile.readline()
		line = line.rstrip('\n')
	infile.close()
	outfile.close()
	return('Done')

def sort_file(args):
	fname_prefix = args[0]
	parnum = args[1]
	os.system('sort --parallel='+str(parnum)+' -o '+fname_prefix+'_sorted '+fname_prefix)                             # make -uo to remove duplicates
	return('Done')

def unique_amino(args):
	fpath = args[0]
	fname_prefix1 = args[1]
	fname_prefix2 = args[2]
	os.system('comm -23 '+fpath+'/'+fname_prefix1+' '+fpath+'/'+fname_prefix2+' > '+fpath+'/'+fname_prefix1+'_unq')
	return('Done')

def expr_cutoff(args):
	fname = args[0]
	cut = args[1]
	os.system('awk \'{if($1>='+cut+'){print i++ " " $2}}\' '+fname+' > '+fname+'_hiexpseqs')
	os.system('awk \'{if($1>='+cut+'){print i++ " " $1}}\' '+fname+' > '+fname+'_hiexpcnts')
	return('Done')

def get_binding(args):
	fname = args[0]
	netMHCpan_path = args[1]
	HLAalleles = args[2]
	fin = open(fname,'r')
	fout = open(fname+'.fa','w')
	line = fin.readline().rstrip('\n')
	while line != '':
		line = line.split(' ')
		outline = line[1]
		ind = outline.find('_')
		if ind > 6:
			fout.write('>amino_'+line[0]+'\n')
			fout.write(outline[0:ind]+'\n')
		elif ind == -1:
			fout.write('>amino_'+line[0]+'\n')
			fout.write(outline+'\n')
		else:
			pass
		line = fin.readline().rstrip('\n')
	fout.close()
	fin.close()
	command_string = 'cat';
	for HLAallele in HLAalleles:
		tmp_dir = fname.split('/')
		tmp_dir = '/tmp_'+tmp_dir[len(tmp_dir)-1]
		os.system(netMHCpan_path+'/netMHCpan_v2 '+tmp_dir+' -a '+HLAallele+' '+fname+'.fa > '+fname+'_affinity_'+HLAallele)
		os.system('rm -rf '+netMHCpan_path+tmp_dir)
		command_string = command_string + ' ' + fname + '_affinity_' + HLAallele
	command_string = command_string + ' > ' + fname+'_affinity'
	os.system(command_string)
	#for HLAallele in HLAalleles:
	#	os.system('rm ' + fname+'_affinity_'+HLAallele)
	return('Done')

def is_number(s):
	try:
		int(s)
		return(True)
	except ValueError:
		return(False)

def clean_binding(args):
	fname = args
	print(fname)
	fin = open(fname,'r')
	fout = open(fname+'_clean','w')
	line = fin.readline();
	while line != '':
		line = line.strip('\n')
		if line == '':
			pass
		elif line[0] == '#':
			pass
		elif line[0] =='-':
			line = fin.readline()
			line = fin.readline()
		elif line[0] == ' ' or line[0] == '\t':
			line = line.split()
			if is_number(line[0]):
				fout.write(",".join(line)+'\n')
			else:
				pass;
		else:
			pass
		line = fin.readline()
	return('Done')

print('Initialization Done: '+str(round(time.time() - start_time,4))+'s')

def main():
	path = '/N/dc2/scratch/johnstrs/TCGA/ControlledData/'
	meta_data = 'gdc_sample_sheet.2019-04-12_LUAD.tsv'
	cancer_type = 'LUAD'
	netMHCpan_path = '/N/u/johnstrs/BigRed2/neoantigen_proj/netMHCpan-4.0'
	convert_bam_to_fastq = True
	gen_uuid_mapping = False
	patient2process = 1
	nfiles = 32
	if gen_uuid_mapping:
		uuid_map(meta_data,cancer_type)
	with open(cancer_type+'_uuids.txt','r') as f:
		i = 1
		line = f.readline().rstrip('\n')
		while i < patient2process:
			line = f.readline().rstrip('\n')
			i = i+1
		line = line.split('\t')
	print( line )
	tumor_fpath = path+cancer_type+'/'+line[1]+'/'+line[2]
	normal_fpath = path+cancer_type+'/'+line[3]+'/'+line[4]
	patseq_path = path+cancer_type+'/'+line[0]
	print(tumor_fpath)
	print(normal_fpath)
	print(patseq_path)
	# convert bam to fastq if necessary
	'''
	if convert_bam_to_fastq:
		bam2fastq((path+cancer_type+'/'+line[3],line[4].split('.',1)[0]))
	# Getting the HLA alleles based on the RNAseq reads
	'''
	hlaPreds = fastq2hla((path+cancer_type+'/'+line[3],line[4].split('.',1)[0],'/N/u/johnstrs/BigRed2/tools',patient2process))
	# Retrieving sequences from files
	pool = mp.Pool(processes=2)
	if not os.path.isdir(patseq_path):
		print('Creating new output directory')
		os.system('mkdir '+patseq_path)
		stats = pool.map(bam2seq,((tumor_fpath,patseq_path+'/tumor.txt'),(normal_fpath,patseq_path+'/normal.txt')))
	else:
		print('Output directory already exists. Is this a duplicate process?')
		print('Deleting old output directory and replacing with new one: '+patseq_path)
		os.system('rm '+patseq_path+'/*')
		stats_bam2seq = pool.map(bam2seq,((tumor_fpath,patseq_path+'/tumor.txt'),(normal_fpath,patseq_path+'/normal.txt')))
	pool = mp.Pool(processes=2)
	stats_split_file1 = pool.map(split_file, ((patseq_path+'/tumor',nfiles/2),(patseq_path+'/normal',nfiles/2)))
	# Begin multiprocessing the subfiles into amino sequences
	pool = mp.Pool(processes=32)
	stats_na2aa = pool.map(translate_file,(patseq_path+'/tumoraa',patseq_path+'/tumorab',patseq_path+'/tumorac',patseq_path+'/tumorad',
	patseq_path+'/tumorae',patseq_path+'/tumoraf',patseq_path+'/tumorag',patseq_path+'/tumorah',
	patseq_path+'/tumorai',patseq_path+'/tumoraj',patseq_path+'/tumorak',patseq_path+'/tumoral',
	patseq_path+'/tumoram',patseq_path+'/tumoran',patseq_path+'/tumorao',patseq_path+'/tumorap',
	patseq_path+'/normalaa',patseq_path+'/normalab',patseq_path+'/normalac',patseq_path+'/normalad',
	patseq_path+'/normalae',patseq_path+'/normalaf',patseq_path+'/normalag',patseq_path+'/normalah',
	patseq_path+'/normalai',patseq_path+'/normalaj',patseq_path+'/normalak',patseq_path+'/normalal',
	patseq_path+'/normalam',patseq_path+'/normalan',patseq_path+'/normalao',patseq_path+'/normalap'))
	pool = mp.Pool(processes=2)
	stats_merge_file = pool.map(merge_files, ((patseq_path,'normal',16,'_amino'),(patseq_path,'tumor',16,'_amino')))
	split_file((patseq_path+'/tumor_amino',32))
	sort_file((patseq_path+'/normal_amino',32))
	pool = mp.Pool(processes=32)
	stats_sort = pool.map(sort_file,((patseq_path+'/tumor_aminoaa',1),(patseq_path+'/tumor_aminoab',1),(patseq_path+'/tumor_aminoac',1),(patseq_path+'/tumor_aminoad',1),
	(patseq_path+'/tumor_aminoae',1),(patseq_path+'/tumor_aminoaf',1),(patseq_path+'/tumor_aminoag',1),(patseq_path+'/tumor_aminoah',1),
	(patseq_path+'/tumor_aminoai',1),(patseq_path+'/tumor_aminoaj',1),(patseq_path+'/tumor_aminoak',1),(patseq_path+'/tumor_aminoal',1),
	(patseq_path+'/tumor_aminoam',1),(patseq_path+'/tumor_aminoan',1),(patseq_path+'/tumor_aminoao',1),(patseq_path+'/tumor_aminoap',1),
	(patseq_path+'/tumor_aminoaq',1),(patseq_path+'/tumor_aminoar',1),(patseq_path+'/tumor_aminoas',1),(patseq_path+'/tumor_aminoat',1),
	(patseq_path+'/tumor_aminoau',1),(patseq_path+'/tumor_aminoav',1),(patseq_path+'/tumor_aminoaw',1),(patseq_path+'/tumor_aminoax',1),
	(patseq_path+'/tumor_aminoay',1),(patseq_path+'/tumor_aminoaz',1),(patseq_path+'/tumor_aminoba',1),(patseq_path+'/tumor_aminobb',1),
	(patseq_path+'/tumor_aminobc',1),(patseq_path+'/tumor_aminobd',1),(patseq_path+'/tumor_aminobe',1),(patseq_path+'/tumor_aminobf',1)))
	# generate unique amino files
	pool = mp.Pool(processes=32)
	stats_tumor_unique = pool.map(unique_amino,((patseq_path,'tumor_aminoaa_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoab_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoac_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoad_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoae_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoaf_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoag_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoah_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoai_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoaj_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoak_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoal_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoam_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoan_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoao_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoap_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoaq_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoar_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoas_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoat_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoau_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoav_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoaw_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoax_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminoay_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoaz_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminoba_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminobb_sorted','normal_amino_sorted'),
	(patseq_path,'tumor_aminobc_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminobd_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminobe_sorted','normal_amino_sorted'),(patseq_path,'tumor_aminobf_sorted','normal_amino_sorted')))
	merge_files((patseq_path,'tumor_amino',32,'_sorted_unq'))
	sort_file((patseq_path+'/tumor_amino_sorted_unq',32))
	os.system('uniq -cd '+patseq_path+'/tumor_amino_sorted_unq_sorted > '+patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates')
	# Generate binding affinity files
	expr_cutoff((patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates',10))
	split_file((patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqs',32))
	print('Predicting Binding Affinities')
	pool = mp.Pool(processes=32)
	stats_affinity = pool.map(get_binding,((patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaa',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsab',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsac',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsad',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsae',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaf',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsag',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsah',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsai',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaj',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsak',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsal',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsam',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsan',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsao',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsap',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaq',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsar',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsas',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsat',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsau',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsav',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaw',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsax',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsay',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsaz',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsba',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsbb',netMHCpan_path,hlaPreds),
	(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsbc',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsbd',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsbe',netMHCpan_path,hlaPreds),(patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqsbf',netMHCpan_path,hlaPreds)))
	print('Done')
	merge_files((patseq_path,'tumor_amino_sorted_unq_sorted_duplicates_hiexpseqs',32,'_affinity'))
	clean_binding((patseq_path+'/tumor_amino_sorted_unq_sorted_duplicates_hiexpseqs_affinity'))
	# Code works up to here
	
	
if __name__ == "__main__":
	main()
	print('Finished: '+str(round(time.time()-start_time,6))+'s')
