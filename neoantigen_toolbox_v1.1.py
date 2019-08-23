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
	print('awk \'{print $3 " " $4 " " $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
	os.system('awk \'{print $3 " " $4 " " $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
	print('rm '+fname[0:fname.find('.')]+'.sam')
	#os.system('rm '+fname[0:fname.find('.')]+'.sam')
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
	print('Splitting '+fname+' into '+str(nfiles)+' files')
	# Begin by counting the lines of input sequence file
	#print('Counting lines')
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
	#print('Done: '+str(round(time.time() - start_time,4))+'s \tLines: '+str(N))
	# Split input sequence file into equal sized subfiles for multiprocessing
	#print('Splitting into multiple files')
	if os.path.isfile(fname+'.txt'):
		os.system('split -l '+str(N/nfiles+1)+' '+fname+'.txt '+fname)
	else:
		os.system('split -l '+str(N/nfiles+1)+' '+fname+' '+fname)
	print('Done splitting '+fname+' into '+str(nfiles)+' files:'+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
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
	print('Done merging ' + fname_prefix + ' files into one file:'+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
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

def write_kmers(line,meta,outfile,k):
    idx = 0;
    ll = len(line);
    while idx+k-1<ll:
        outfile.write(str(meta[0])+' '+str(meta[1])+' '+str(meta[2])+' '+line[idx:idx+k]+'\n');
        idx = idx+1;
    return('Done')

def translate_file(args):
    #print('Translating '+fname_prefix_sub+' to '+str(k)+'mers ')
    fname_prefix_sub = args[0]
    k = args[1]
    infile = open(fname_prefix_sub, 'r')
    outfile = open(fname_prefix_sub+'_aminodat', 'w')
    i=0;
    line = infile.readline()
    line = line.rstrip('\n')
    while line != '':
        line = line.split(' ')
        if 'N' not in line[2]:
            #print('found')
            tmp = line[2][0:(len(line[2])-(len(line[2])%3))]
            #print(tmp)
            aa1 = translate(tmp)
            tmp = line[2][1:(len(line[2])-((len(line[2])-1)%3))]
            aa2 = translate(tmp)
            tmp = line[2][2:(len(line[2])-((len(line[2])-2)%3))]
            aa3 = translate(tmp)
            write_kmers(aa1,line,outfile,k)
            write_kmers(aa2,line,outfile,k)
            write_kmers(aa3,line,outfile,k)
        line = infile.readline()
        line = line.rstrip('\n')
    infile.close()
    outfile.close()
    return('Done')

def get_col(args):
	fpath = args[0]
	fnameprefix = args[1]
	fnamepostfix = args[2]
	onamepostfix = args[3]
	col = args[4]
	os.system('awk \'{print $'+str(col)+'}\' '+fpath+'/'+fnameprefix+fnamepostfix+' > '+fpath+'/'+fnameprefix+onamepostfix)
	return('Done')

def sort_file(args):
	fpath = args[0]
	fname = args[1]
	parnum = args[2]
	os.system('sort -f -T '+fpath+' -d --parallel='+str(parnum)+' -o '+fpath+'/'+fname+'_sorted '+fpath+'/'+fname)                             # make -uo to remove duplicates
	return('Done')

def unique_amino(args):
	fpath = args[0]
	fname_prefix1 = args[1]
	fname_prefix2 = args[2]
	os.system('comm --check-order -23 '+fpath+'/'+fname_prefix1+' '+fpath+'/'+fname_prefix2+' > '+fpath+'/'+fname_prefix1+'_unq')
	return('Done')

def expr_cutoff(args):
	fname = args[0]
	cut = args[1]
	cmd = 'awk \' $1 >= '+str(cut)+' \' '+fname+' > '+fname+'_hiexp'
	cmd1 = 'awk \'{if($1>'+str(cut)+'){print $2}}\' '+fname+' > '+fname+'_hiexpseqs'
	#cmd2 = 'awk \'{if($1>='+str(cut)+'){print i++ " " $1}}\' '+fname+' > '+fname+'_hiexpcnts'
	print(cmd)
	print(cmd1)
	#print(cmd2)
	os.system(cmd)
	os.system(cmd1)
	#os.system(cmd2)
	return('Done')

def get_binding(args):
	fname = args[0]
	netMHCpan_path = args[1]
	HLAalleles = args[2]
	kmer_k = args[3]
	'''
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
	'''
	command_string = 'cat';
	for HLAallele in HLAalleles:
		tmp_dir = fname.split('/')
		tmp_dir = '/tmp_'+tmp_dir[len(tmp_dir)-1]
		os.system(netMHCpan_path+'/netMHCpan_v2 '+tmp_dir+' -a '+HLAallele+' -l '+str(kmer_k)+' -p '+fname+' > '+fname+'_affinity_'+HLAallele)
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

def subset_metadata(args):
    dbfname = args[0];
    queryfname = args[1];
    dbinfile = open(dbfname,'r')
    query_seqs = [line.rstrip('\n') for line in open(queryfname,'r')]
    print(dbfname)
    print(query_seqs)
    outfile = open(dbfname+'_metadata','w')
    line = dbinfile.readline()
    line2 = line.rstrip('\n')
    while line != '':
        line2 = line2.split(' ')
        if line2[3] in query_seqs:
            outfile.write(line)
        line = dbinfile.readline()
        line2 = line.rstrip('\n')
    return('Done')

print('Initialization Done: '+str(round(time.time() - start_time,4))+'s')

def main():
	path = '/N/dc2/scratch/johnstrs/TCGA/ControlledData/'
	meta_data = 'gdc_sample_sheet.2019-04-12_LUAD.tsv'
	cancer_type = 'LUAD'
	netMHCpan_path = '/N/u/johnstrs/BigRed2/neoantigen_proj/netMHCpan-4.0'
	convert_bam_to_fastq = True
	gen_uuid_mapping = False
	patient2process = int(sys.argv[1])
	nfiles = 32
	kmer_k = int(sys.argv[2])
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
	if convert_bam_to_fastq:
		bam2fastq((path+cancer_type+'/'+line[3],line[4].split('.',1)[0]))
	# Getting the HLA alleles based on the RNAseq reads
	print('GETTING HLA ALLELES')
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
	print('TRANSLATING FILES TO AMINO ACIDS')
	pool = mp.Pool(processes=2)
	stats_split_file1 = pool.map(split_file, ((patseq_path+'/tumor',nfiles/2),(patseq_path+'/normal',nfiles/2)))
	# Begin multiprocessing the subfiles into amino sequences
	pool = mp.Pool(processes=32)
	stats_na2aa = pool.map(translate_file,((patseq_path+'/tumoraa',kmer_k),(patseq_path+'/tumorab',kmer_k),(patseq_path+'/tumorac',kmer_k),(patseq_path+'/tumorad',kmer_k),
	(patseq_path+'/tumorae',kmer_k),(patseq_path+'/tumoraf',kmer_k),(patseq_path+'/tumorag',kmer_k),(patseq_path+'/tumorah',kmer_k),
	(patseq_path+'/tumorai',kmer_k),(patseq_path+'/tumoraj',kmer_k),(patseq_path+'/tumorak',kmer_k),(patseq_path+'/tumoral',kmer_k),
	(patseq_path+'/tumoram',kmer_k),(patseq_path+'/tumoran',kmer_k),(patseq_path+'/tumorao',kmer_k),(patseq_path+'/tumorap',kmer_k),
	(patseq_path+'/normalaa',kmer_k),(patseq_path+'/normalab',kmer_k),(patseq_path+'/normalac',kmer_k),(patseq_path+'/normalad',kmer_k),
	(patseq_path+'/normalae',kmer_k),(patseq_path+'/normalaf',kmer_k),(patseq_path+'/normalag',kmer_k),(patseq_path+'/normalah',kmer_k),
	(patseq_path+'/normalai',kmer_k),(patseq_path+'/normalaj',kmer_k),(patseq_path+'/normalak',kmer_k),(patseq_path+'/normalal',kmer_k),
	(patseq_path+'/normalam',kmer_k),(patseq_path+'/normalan',kmer_k),(patseq_path+'/normalao',kmer_k),(patseq_path+'/normalap',kmer_k)))
	print('RETRIEVING AMINO ACIDS')
	pool = mp.Pool(processes=32)
	stats_getamino = pool.map(get_col,((patseq_path,'tumoraa','_aminodat','_amino',4),(patseq_path,'tumorab','_aminodat','_amino',4),(patseq_path,'tumorac','_aminodat','_amino',4),(patseq_path,'tumorad','_aminodat','_amino',4),
	(patseq_path,'tumorae','_aminodat','_amino',4),(patseq_path,'tumoraf','_aminodat','_amino',4),(patseq_path,'tumorag','_aminodat','_amino',4),(patseq_path,'tumorah','_aminodat','_amino',4),
	(patseq_path,'tumorai','_aminodat','_amino',4),(patseq_path,'tumoraj','_aminodat','_amino',4),(patseq_path,'tumorak','_aminodat','_amino',4),(patseq_path,'tumoral','_aminodat','_amino',4),
	(patseq_path,'tumoram','_aminodat','_amino',4),(patseq_path,'tumoran','_aminodat','_amino',4),(patseq_path,'tumorao','_aminodat','_amino',4),(patseq_path,'tumorap','_aminodat','_amino',4),
	(patseq_path,'normalaa','_aminodat','_amino',4),(patseq_path,'normalab','_aminodat','_amino',4),(patseq_path,'normalac','_aminodat','_amino',4),(patseq_path,'normalad','_aminodat','_amino',4),
	(patseq_path,'normalae','_aminodat','_amino',4),(patseq_path,'normalaf','_aminodat','_amino',4),(patseq_path,'normalag','_aminodat','_amino',4),(patseq_path,'normalah','_aminodat','_amino',4),
	(patseq_path,'normalai','_aminodat','_amino',4),(patseq_path,'normalaj','_aminodat','_amino',4),(patseq_path,'normalak','_aminodat','_amino',4),(patseq_path,'normalal','_aminodat','_amino',4),
	(patseq_path,'normalam','_aminodat','_amino',4),(patseq_path,'normalan','_aminodat','_amino',4),(patseq_path,'normalao','_aminodat','_amino',4),(patseq_path,'normalap','_aminodat','_amino',4)))
	print('MERGING TUMOR AND NORMAL AMINO FILES')
	pool = mp.Pool(processes=4)
	stats_merge_file = pool.map(merge_files, ((patseq_path,'normal',16,'_aminodat'),(patseq_path,'tumor',16,'_aminodat'),(patseq_path,'normal',16,'_amino'),(patseq_path,'tumor',16,'_amino')))
	print('SORTING TUMOR AND NORMAL AMINO FILES')
	sort_file((patseq_path,'tumor_amino',8))
	sort_file((patseq_path,'normal_amino',8))
	print('GETTING UNIQUE SEQUENCES FOR TUMOR AND NORMAL AMINO')
	os.system('uniq -cd '+patseq_path+'/tumor_amino_sorted > '+patseq_path+'/tumor_amino_sorted_dups')
	os.system('uniq -cd '+patseq_path+'/normal_amino_sorted > '+patseq_path+'/normal_amino_sorted_dups')
	print('REMOVING SEQUENCES BELOW EXPRESSION CUTOFF')
	expr_cutoff((patseq_path+'/tumor_amino_sorted_dups',20))
	expr_cutoff((patseq_path+'/normal_amino_sorted_dups',5))
	print('SPLITTING TUMOR AMINO FILES')
	split_file((patseq_path+'/tumor_amino_sorted_dups_hiexpseqs',32))
	# generate unique amino files
	print('GETTING UNIQUE AMINO ACID SEQUENCES')
	pool = mp.Pool(processes=32)
	stats_tumor_unique = pool.map(unique_amino,((patseq_path,'tumor_amino_sorted_dups_hiexpseqsaa','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsab','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsac','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsad','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsae','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsaf','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsag','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsah','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsai','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsaj','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsak','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsal','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsam','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsan','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsao','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsap','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsaq','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsar','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsas','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsat','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsau','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsav','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsaw','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsax','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsay','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsaz','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsba','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsbb','normal_amino_sorted_dups_hiexpseqs'),
	(patseq_path,'tumor_amino_sorted_dups_hiexpseqsbc','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsbd','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsbe','normal_amino_sorted_dups_hiexpseqs'),(patseq_path,'tumor_amino_sorted_dups_hiexpseqsbf','normal_amino_sorted_dups_hiexpseqs')))
	merge_files((patseq_path,'tumor_amino_sorted_dups_hiexpseqs',32,'_unq'))
	print('PREDICTING BINDING AFFINITIES')
	stats_affinity = get_binding((patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq',netMHCpan_path,hlaPreds,kmer_k))
	clean_binding((patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq_affinity'))
	print('GETTING METADATA FOR UNIQUE KMERS')
	pool = mp.Pool(processes=16)
	stats_metadata = pool.map(subset_metadata,((patseq_path+'/tumoraa_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorab_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorac_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorad_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	(patseq_path+'/tumorae_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoraf_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorag_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorah_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	(patseq_path+'/tumorai_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoraj_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorak_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoral_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	(patseq_path+'/tumoram_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoran_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorao_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorap_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq')))
	#(patseq_path+'/tumoraq_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorar_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoras_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorat_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	#(patseq_path+'/tumorau_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorav_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoraw_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorax_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	#(patseq_path+'/tumoray_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumoraz_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorba_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorbb_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),
	#(patseq_path+'/tumorbc_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorbd_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorbe_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq'),(patseq_path+'/tumorbf_aminodat',patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq')))
	merge_files((patseq_path,'tumor',16,'_aminodat_metadata'))
	print('CLEANING DIRECTORY AND WRITING FILES')
	os.system('mv '+patseq_path+'/tumor_aminodat_metadata '+path+cancer_type+'/'+line[0]+'_'+str(patient2process)+'_'+str(kmer_k)+'mer_metadata.txt')
	os.system('mv '+patseq_path+'/tumor_amino_sorted_dups_hiexpseqs_unq_affinity_clean '+path+cancer_type+'/'+line[0]+'_'+str(patient2process)+'_'+str(kmer_k)+'mer_neoantigens.txt')
	os.system('mv '+patseq_path+'/tumor_amino_sorted_dups '+path+cancer_type+'/'+line[0]+'_'+str(patient2process)+'_'+str(kmer_k)+'mer_tumorpep.txt')
	os.system('mv '+patseq_path+'/normal_amino_sorted_dups '+path+cancer_type+'/'+line[0]+'_'+str(patient2process)+'_'+str(kmer_k)+'mer_normalpep.txt')
	os.system('rm -rf '+patseq_path)
	

if __name__ == "__main__":
	main()
	print('Finished: '+str(round(time.time()-start_time,6))+'s')
