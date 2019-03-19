#! /usr/bin/python

# Copyright (c) Sara Vieira-Silva & Eduardo Rocha 
#				Institut Pasteur - Paris - France
#				Atelier de BioInformatique Universite Pierre et Marie Curie - Paris - France
#
# desc: calculates minimum generation time for prokaryotes based on codon usage bias intensity
#
# History: Oct 2009: 1.00: single organism, bootstrap on utilCUB results
# History: Dec 2009: 1.01: single organism with bootstrap on sequences or metagenome with bootstrap on utilCUB results, blast available to retrieve HEG
# History: Dec 2009: 1.02: deprecated os.system() replaced with subprocess()
#History: Dec 2009: 1.03: example files option added
#History: Dec 2009: 1.04: empty sequences bug removed (with identifier bu len=0), test valid OGT
#History: Dec 2009: 1.05: getOGT debugged, test for valid cub data or error in R script
#History: Jan 2010: 1.06: getOGT debugged, comma/dots in awk
#history:  Jan 2010: 1.07: example outputs updated

GROWTHPRED_VERSION = "GROWTHPRED VERSION 1.07"

import math
import sys
import string
import os
import subprocess
import random
from optparse import OptionParser

GROWTHPRED_DIR = os.path.dirname(os.path.abspath(__file__))
GROWTHPRED_BIN = GROWTHPRED_DIR + "/"
GROWTHPRED_LIBEX = os.path.join(GROWTHPRED_DIR, "Programs/")
GROWTHPRED_SHARE = os.path.join(GROWTHPRED_DIR, "shared/")
GROWTHPRED_TMP = "./"
GROWTHPRED_R = GROWTHPRED_BIN

'''
try:
	GROWTHPRED_BIN = os.environ['GROWTHPRED_BIN'] #/local/gensoft/bin
except:
	GROWTHPRED_BIN = ""
if GROWTHPRED_BIN not in sys.path:
	sys.path.append( GROWTHPRED_BIN )

try:	
	GROWTHPRED_LIBEX = os.environ['GROWTHPRED_LIBEX'] #/local/gensoft/libexec/growthpred
except:
	GROWTHPRED_LIBEX = "./Programs/"
if GROWTHPRED_LIBEX not in sys.path:
	sys.path.append( GROWTHPRED_LIBEX )
	
try:
	GROWTHPRED_SHARE = os.environ['GROWTHPRED_SHARE'] #/local/gensoft/share/program
except:
	GROWTHPRED_SHARE = "./shared/" 
if GROWTHPRED_SHARE not in sys.path:
	sys.path.append( GROWTHPRED_SHARE )
	
try:
	GROWTHPRED_TMP = os.environ['GROWTHPRED_TMP'] #directory created for each job
except:
	GROWTHPRED_TMP = "./" 
if GROWTHPRED_TMP not in sys.path:
	sys.path.append( GROWTHPRED_TMP )

try:
       GROWTHPRED_R = os.environ['GROWTHPRED_R']
except:
       GROWTHPRED_R = GROWTHPRED_BIN
'''
		
def Parsopt():
	print GROWTHPRED_VERSION
	if (len(sys.argv[0:]) < 2):
		print "type \"growthpred.py -h\ for help"
		os._exit(0)
	parser = OptionParser()
	parser.add_option("-f","--heg", dest="HEG", help="File of highly expressed genes (DNA, FASTA)")
	parser.add_option("-g","--nonheg", dest="OTH", help="Complete genome (DNA, FASTA)")
        parser.add_option("-d","--dirinput", dest="dirinput", help="input directory where heg and nonheg can be found")
	parser.add_option("-o", dest="outfile", help="Name of output file")
	parser.add_option("-c", action="store", dest="code", type="int", help="0 : Universal; 1 : Yeast mitochondrial; 2 : Vertebrate mitochondrial; 3 : Mold/Protozoan/Mycoplasma/Spiroplasma; 4 : Invertebrate mitochondrial; 5 : Candida cylindracea; 6 : Ciliate; 7 : Euplotes; 8 : Echinoderm mitochondrial")
	parser.add_option("-s", dest="start", action="store_true", default=False, help="Ignore first (start) codons")
	parser.add_option("-S", dest="stop", action="store_true", default=False, help="Ignore last  (stop)  codons")
	parser.add_option("-T", action="store", dest="OGT", type="float", default=36, help="Optimal growth temperature (Celsius)")
	parser.add_option("-t", dest="tpred", action="store_true", default=False, help="Estimate optimal growth temperature")
	parser.add_option("-b", dest="blast", action="store_true", default=False, help="Retrieve ribosomal protein genes by blast")
	parser.add_option("-r", dest="filerp", action="store_true", default=False, help="Recover file with ribosomal protein genes retrieved by blast")
	parser.add_option("-i", dest="filecub", action="store_true", default=False, help="Recover file with codon usage bias indexes for each gene")
	parser.add_option("-m", dest="metagen", action="store_true", default=False, help="Metagenome or mixed organisms sequences")
	parser.add_option("-e", dest="example", action="store_true", default=False, help="Use E. coli K12 example files")
	return parser.parse_args()

def Liststops(code):
	#dictionnary of stop codons for each genetic code
	dictstops = {
				'0':['TAA','TAG','TGA'], \
				'1':['TAA','TAG'], \
				'2':['TAA','TAG','AGA','AGG'], \
				'3':['TAA','TAG'], \
				'4':['TAA','TAG'], \
				'5':['CTA','TAA','TAG','TGA'], \
				'6':['TGA'], \
				'7':['TAA','TAG'], \
				'8':['TAA','TAG'] \
				}
	
	stopcods = dictstops[code]
	return stopcods
	
def CheckStopCodons(dnaSeq, code):
	#dnaSeq length is multiple of 3, no header, no eols
	stopCodons = Liststops(code)
	if ("\n" in dnaSeq):
		dnaSeq = dnaSeq.replace("\n","")
	n = len(dnaSeq)
	stops_exist = False
	for i in xrange(0,n,3):
		if (dnaSeq[i:i+3] in stopCodons):
			stops_exist = True
			break
	
	return stops_exist 

def Getseqs(dirdata, file, dirtmp, (options, args)):
	#input sequences file is turned into dictionnary, removes starts and stops if opted, filters out seqs non mltiple of 3 or with internal stops
	dico = {}	
	f1 = open("%s%s" %(dirdata,file),"r")
	lines = f1.readlines()
	f1.close()
	numseq = 0
	seq = ""
	n=0
	fileEr = open("%s%s.errors" %(dirtmp,options.outfile),"a")
	while n < len(lines):
		if (n < len(lines)) and (lines[n][0] == ">"):
			seq = ""
			pt = string.split(lines[n])[0]
			k = n + 1
			numseq = numseq + 1
			while (k < len(lines)) and (lines[k][0] != ">"):
				seq = seq + lines[k].rstrip()
				k = k + 1
			n = k
			name = pt[1:]
			if (len(seq)%3 == 0):
				if (options.stop == True) & (options.start == False):
					seq=seq[:-3]
				elif (options.stop == True) & (options.start == True):
					seq=seq[3:-3]
				elif (options.stop == False) & (options.start == True):
					seq=seq[3:]
				else:
					seq=seq
				if (CheckStopCodons(seq, str(options.code)) == False):
					if (len(seq) != 0): 
						dico[name]=seq
					else:
						fileEr.write("ERROR: empty sequence: %s \n" %(name))
				else:
					fileEr.write("ERROR: internal stop in sequence: %s \n" %(name))
			else:
				fileEr.write("ERROR: sequence length not multiple of 3: %s \n" %(name))
		else:
			n = n + 1	
	fileEr.close()
	return (dico, numseq) 
	
def Bootseqs(dirfile, file, niteration, dico,type):
	# x concatenation of n randomly picked sequences  (x = niteration; n = number of sequences in file)
	print "Bootstrap", type
	tempout = open("%s%s" %(dirfile, file),"a")
	for randi in range(0,niteration):
		seq = ""
		name = "seq"+file+str(randi)
		for nseqs in range(0,len(dico)):
			seq = seq+random.choice(dico.values())
		tempout.write(">%s \n%s \n" %(name,seq))
	tempout.close()
	
def Getnucfreq(direx, dirfile, file):
	#calculates nucleotide frequency in quartets (neutral selection)
	nucfreqs = subprocess.Popen(["%swcountq" %(direx), "%s%s" %(dirfile,file)], stdout=subprocess.PIPE).communicate()[0]
	A = str.split(nucfreqs)[45]
 	C = str.split(nucfreqs)[48]
 	G = str.split(nucfreqs)[51]
	return (A,C,G)
			
def GetCUB(direx,A,C,G,code,dirfile,file,type,outfile):
	#calculates codon usage bias indexes ENCp and P (Vieira-Silva and Rocha, PlosGen)
	print "Calculate codon usage bias indexes for %s sequences" %(type)
        cmd='%sutilCUB -A %s -C %s -G %s -c %s -t < %s%s' %(direx,A,C,G,code,dirfile,file)
        print(cmd)
	p1 = subprocess.Popen('%sutilCUB -A %s -C %s -G %s -c %s -t < %s%s' %(direx,A,C,G,code,dirfile,file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	f1 = open("%s%s.errors" %(dirfile,outfile),"a")
	f1.write(p1[1])
	f1.close()  
	f2 = open("%s%s.%s.enc.si" %(dirfile,outfile,type),"w")
	f2.write(p1[0])
	f2.close() 
	
def Getgentime_single(direx,dirshare,dirfile,outfile,OGT):
	#calculates minimum generation time (Vieira-Silva and Rocha, PlosGen)
	print "Calculate predicted minimum generation time"
	p1 = subprocess.Popen('R --slave -q --no-save < %sgetd_single.R %s%s.OTH.enc.si %s%s.HEG.enc.si %s %s%s.results' %(dirshare,dirfile,outfile,dirfile,outfile,OGT,dirfile,outfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	f1 = open("%s%s.errors" %(dirfile,outfile),"a")
	f1.write(p1[1])
	f1.close()  
	
def Getgentime_mixed(direx,dirshare,dirfile,outfile,OGT):
	#calculates minimum generation time (Vieira-Silva and Rocha, PlosGen)
	print "Calculate predicted minimum generation time"
	p1 = subprocess.Popen('R --slave -q --no-save < %sgetd_mixed.R %s%s.OTH.enc.si %s%s.HEG.enc.si %s %s%s.results' %(dirshare,dirfile,outfile,dirfile,outfile,OGT,dirfile,outfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	f1 = open("%s%s.errors" %(dirfile,outfile),"a")
	f1.write(p1[1])
	f1.close() 	
	
def Getogt(direx, dirfile, file, code, dirbin):
	#calculates optimal growth temperature (Zeldovich et al, PlosCompBiol, 2007)
	print "Calculate predicted optimal growth temperature"
	ogt = subprocess.Popen(["%sget_OGT.sh" %(direx), "%s" %(code), "%s%s" %(dirfile,file),"%s" %(dirbin)], stdout=subprocess.PIPE, stderr=None).communicate()[0]
	return ogt

def Getheg(direx, dirfile, file, hegdb,diroutfile,outfile,blastall,getentryF):
	#retrieves the set of seq from input file that have high similarity (Evalue<0.00001, identity>40%) and similar length (+-10%) to ribosomal database seqs 
	print "Blast against ribosomal protein database"
	subprocess.Popen('%sget_rprots.sh %s%s %s %s%s %s %s' %(direx,dirfile,file,hegdb,diroutfile,outfile,blastall,getentryF), shell=True, stdout=None, stderr=None).communicate()[0]
	
def Getcai(diremb,dirfile,fileheg,fileoth):
	#calculates codon adaptation index (Sharp and LI, NAR, 1987)
	subprocess.Popen('%scusp %s%s %s%s.cusp' %(diremb,dirfile,fileheg,dirfile,fileheg), shell=True, stdout=None, stderr=None).communicate()
	subprocess.Popen('%scai %s%s -cfile %s%s.cusp %s%s.cai' %(diremb,dirfile,fileoth,dirfile,fileheg,dirfile,fileoth), shell=True, stdout=None, stderr=None).communicate()
	subprocess.Popen('rm -f %s%s.cusp' %(dirfile,fileheg), shell=True, stdout=None, stderr=None).communicate()
	
def Join(dirfile,fileenc,filecai,outfile):
	#Creates CUB file by joining ENCp+P file with CAI file
	print "Create codon usage bias indexes file"
	p1=subprocess.Popen('cut -f 2,4 -d " " %s%s' %(dirfile,filecai), shell=True, stdout=subprocess.PIPE)
	p2=subprocess.Popen('sort', shell=True, stdout=subprocess.PIPE, stdin=p1.stdout).communicate()
	f1 = open("%s%s.2" %(dirfile,filecai),"w")
	f1.write(p2[0])
	f1.close()  
	p1=subprocess.Popen('sort %s%s' %(dirfile,fileenc), shell=True, stdout=subprocess.PIPE).communicate()
	f1 = open("%s%s.2" %(dirfile,fileenc),"w")
	f1.write(p1[0])
	f1.close()
	f1 = open("%s%s" %(dirfile,outfile),"a")
	f1.write("# NB: The Codon Adaptation Index is only an appropriate proxy \n# for the expression level of genes in fast growing cells \n# minimal doubling times under ~2h \n")
	f1.write("SeqID ENC ENCp P CAI \n")
	p1=subprocess.Popen('join %s%s %s%s' %(dirfile,fileenc+".2",dirfile,filecai+".2"), shell=True, stdout=subprocess.PIPE).communicate()
	f1.write(p1[0])  
	f1.close()
	subprocess.Popen('rm -f %s%s %s%s %s%s %s%s' %(dirfile,fileenc,dirfile,filecai,dirfile,fileenc+".2",dirfile,filecai+".2"), shell=True, stdout=None, stderr=None).communicate()[0]
	
def main():
	
### DEFINE DIRECTORIES
	dirbin = GROWTHPRED_LIBEX #location for utilCUB (C), wcountq (C), translate (C), util_tab_aa (C), getentryF (C)
	dirR = GROWTHPRED_R #location for R
	diremboss = GROWTHPRED_BIN #location for cusp and cai from EMBOSS
	#dirblast = GROWTHPRED_BIN # location for blastall
	dirshare = GROWTHPRED_LIBEX # location for  get_OGT.sh (tcsh), get_rprots.sh (tcsh), getd_mixed.R (R), getd_single.R (R)
	dirdb = GROWTHPRED_SHARE #location for RPROTDB/RPDB/ and example files
	dirtmp = GROWTHPRED_TMP # location for intermediate and output files
	dirdata = GROWTHPRED_TMP # location for input files	
        #print(dirbin, dirR, diremboss, dirblast, dirshare, dirdb, dirtmp, dirdata)

### PARSING OPTIONS	
	(options, args) = Parsopt()	
	if (options.start == True):
		optstart = "-s"
	else:
		optstart = ""
	if (options.stop == True):
		optstop = "-S"
	else:
		optstop = ""
		
	### DATASET - CHOOSE INPUT FILES
	if (options.example == False):
		hegfile = options.HEG
		othfile = options.OTH
                # add dirinput option
		#dirinput = dirdata
                dirinput = options.dirinput + "/"
	else:
		hegfile = "ecoli_ribosomal_genes"
		othfile = "ecoli_complete_genome"
		dirinput = dirdb
	
	### DATASET - CREATE DICTIONARIES
	if (options.blast == True):
		#Getheg(dirshare, dirinput, othfile, dirdb+"RPROTDB/RPDB",dirtmp,options.outfile,dirblast+"blastall",dirbin+"getentryF")
                Getheg(dirshare, dirinput, othfile, dirdb+"RPROTDB/RPDB",dirtmp,options.outfile,"blastall",dirbin+"getentryF")
		getseqsHEG = Getseqs(dirtmp, options.outfile+".ribs", dirtmp, (options, args))
		dicoHEG = getseqsHEG[0]
		numheg = getseqsHEG[1]
		getseqsOTH = Getseqs(dirinput, othfile, dirtmp, (options, args))
		dicoOTH = getseqsOTH[0]
		numoth = getseqsOTH[1]
		if (options.filerp == False):
			subprocess.Popen('rm -f %s%s' %(dirtmp, options.outfile+".ribs"), shell=True, stdout=None, stderr=None).communicate()[0]
	else:
		getseqsHEG = Getseqs(dirinput, hegfile, dirtmp, (options, args))
		dicoHEG = getseqsHEG[0]
		numheg = getseqsHEG[1]
		getseqsOTH = Getseqs(dirinput, othfile, dirtmp, (options, args))
		dicoOTH = getseqsOTH[0]
		numoth = getseqsOTH[1]
	
	if (len(dicoHEG) < 10):
		if (len(dicoHEG) == 0):
			print "No valid highly expressed gene sequence"
			os._exit(0)
		else:
			print "(!) Low number of valid highly expressed gene sequence",len(dicoHEG) 
	if (len(dicoOTH) < 10):
		if (len(dicoOTH) == 0):
			print "No valid control gene sequence"
			os._exit(0)
		else:
			print "(!) Low number of valid control gene sequence:",len(dicoOTH)
	
	###	CALCULATE CODON USAGE BIASES FILE
	if (options.filecub == True) | (options.metagen == True):
		# create clean sequence files for HEG and OTH
		foth = open("%s%s.oths" %(dirtmp,options.outfile),"a")
		fheg = open("%s%s.hegs" %(dirtmp,options.outfile),"a")
		for keyi in dicoHEG.keys():
			fheg.write(">%s\n%s\n" %(keyi,dicoHEG[keyi]))
		fheg.close()
		for keyi in dicoOTH.keys():
			foth.write(">%s\n%s\n" %(keyi,dicoOTH[keyi]))
		foth.close()
		### DATA IS FROM SINGLE ORGANISM	
		if (options.metagen == False):
			# calculate cubs
			nucf=Getnucfreq(dirbin, dirtmp, options.outfile+".oths")
			GetCUB(dirbin,nucf[0],nucf[1],nucf[2],options.code,dirtmp,options.outfile+".oths","oths",options.outfile)
			Getcai(diremboss,dirtmp,options.outfile+".hegs",options.outfile+".oths")
			Join(dirtmp,options.outfile+".oths.enc.si",options.outfile+".oths.cai",options.outfile+".cub")
			subprocess.Popen('rm -f %s%s.hegs %s%s.oths' %(dirtmp,options.outfile,dirtmp,options.outfile), shell=True, stdout=None, stderr=None).communicate()[0]

		### DATA IS FROM METAGENOME	
		else:
			print "Processing metagenome (mixed organisms) sequences"
			if (options.filecub == True):
				print "(!) CAI can only be calculated for single organism datasets"
			nucfreqs=Getnucfreq(dirbin, dirtmp, options.outfile+".oths")
			GetCUB(dirbin,nucfreqs[0],nucfreqs[1],nucfreqs[2],options.code,dirtmp,options.outfile+".oths","OTH",options.outfile)
			GetCUB(dirbin,nucfreqs[0],nucfreqs[1],nucfreqs[2],options.code,dirtmp,options.outfile+".hegs","HEG",options.outfile)
			subprocess.Popen('cp -f %s%s %s%s' %(dirtmp,options.outfile+".OTH.enc.si",dirtmp,options.outfile+".cub"), shell=True, stdout=None, stderr=None).communicate()[0]
                        # test no clean
			subprocess.Popen('rm -f %s%s.hegs %s%s.oths' %(dirtmp,options.outfile,dirtmp,options.outfile), shell=True, stdout=None, stderr=None).communicate()[0]
			
	### DATA IS FROM SINGLE ORGANISM	
	if (options.metagen == False):
		print "Processing single genome sequences"
		
		###	DEFINE N ITERATION BOOTSTRAP
		niteration = 100
		
		### CREATE N RANDOM HEG AND OTH FILES (BOOTSTRAPPED)	
		HEGtmp = options.outfile+".HEG.tmp"
		OTHtmp = options.outfile+".OTH.tmp"
		
		Bootseqs(dirtmp, HEGtmp, niteration, dicoHEG, "HEG sequences")
		Bootseqs(dirtmp, OTHtmp, niteration, dicoOTH, "control sequences")
		
		### CALCULATE QUARTET NUCLEOTIDE FREQUENCY
		nucfreqs=Getnucfreq(dirbin, dirinput, othfile)
		
		### CALCULATE CODON USAGE BIAS
		GetCUB(dirbin,nucfreqs[0],nucfreqs[1],nucfreqs[2],options.code,dirtmp,OTHtmp,"OTH",options.outfile)
		GetCUB(dirbin,nucfreqs[0],nucfreqs[1],nucfreqs[2],options.code,dirtmp,HEGtmp,"HEG",options.outfile)
		
	### CALCULATE OPTIMAL GROWTH TEMPERATURE
	if (options.tpred == True):
		tempC = Getogt(dirshare, dirinput, othfile, options.code,dirbin)
		ogt = round(float(tempC))
		if (options.metagen == True):
			print "Optimal growth temperature prediction is not recommended on mixed organisms sequences / metagenomes"
	else:
		ogt = options.OGT
	if (ogt < 0) | (ogt > 120):
		print "Invalid optimal growth temperature: %sC replaced by default (valid range: 0-120C)" %(ogt)
		ogt = 36
	### CALCULATING	MINIMUM GENERATION TIME	
	if (options.metagen == False):	
		### DATA IS FROM SINGLE ORGANISM
		Getgentime_single(dirR,dirshare,dirtmp,options.outfile,ogt)
	else:
		### DATA IS FROM METAGENOME
		Getgentime_mixed(dirR,dirshare,dirtmp,options.outfile,ogt)
		
	subprocess.Popen('rm -f %s%s.OTH.enc.si %s%s.HEG.enc.si %s%s.OTH.tmp %s%s.HEG.tmp' %(dirtmp,options.outfile,dirtmp,options.outfile,dirtmp,options.outfile,dirtmp,options.outfile), shell=True, stdout=None, stderr=None).communicate()[0]

	
	### FINISH OUTPUT	
	fout = open("%s%s.results" %(dirtmp,options.outfile),"a")
	fout.write("Nucleotide frequencies: A = %0.2f ; C = %0.2f ; G = %0.2f ; T = %0.2f \n" %(float(nucfreqs[0]),float(nucfreqs[1]),float(nucfreqs[2]), 1-(round(float(nucfreqs[0]),2)+round(float(nucfreqs[1]),2)+round(float(nucfreqs[2]),2))))
	fout.write("Input dataset: %s highly expressed genes + %s non-highly expressed genes \n" %(numheg, numoth))
	fout.write("Input optimal growth temperature (Celsius): %s \n" %(ogt))
	fout.close()
	
	### REMOVE REPETITION FROM ERROR
	p1=subprocess.Popen('sort -u %s%s.errors' %(dirtmp,options.outfile), shell=True, stdout=subprocess.PIPE).communicate()
	f1 = open("%s%s.e.2" %(dirtmp,options.outfile),"w")
	f1.write(p1[0])
	subprocess.Popen('mv -f %s%s.e.2 %s%s.errors' %(dirtmp,options.outfile,dirtmp,options.outfile), shell=True, stdout=subprocess.PIPE).communicate()
	f1.close()

main()		


