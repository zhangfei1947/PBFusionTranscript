#!/usr/bin/python
##coding=utf-8
#
# author
#
#   zhangfei@novogene.com
#   2017.02.23
#
# description
#
#   Identify fusion transcript
#
# Usage
#
#   python run_fusion.py -h

import argparse,sys,re,os,glob
from SpliceGrapher.formats.loader import loadGeneModels
from SpliceGrapher.formats.GeneModel import *
from SpliceGrapher.formats import fasta
from SpliceGrapher.SpliceGraph       import *
from SpliceGrapher.shared.GeneModelConverter import *


#------------------
parser = argparse.ArgumentParser(description="Identify fusion transcript")
parser.add_argument('--fasta',help='fasta',required=True)
parser.add_argument('--gmap',help='gmap',required=True)
parser.add_argument('--gmap_idx',help='gmap_idx',required=True)
parser.add_argument('--gmap_idx_p',help='gmap_idx_p',required=True)
parser.add_argument('--gtf',help='gtf',required=True)
parser.add_argument('--STAR',help='STAR')
parser.add_argument('--clean_fq1',help='clean_fq1')
parser.add_argument('--clean_fq2',help='clean_fq2')
#--- to generate karyotype file
parser.add_argument('--genome',help='genome')
parser.add_argument('--abbr',help='abbr')
parser.add_argument('--num',help='num')
#--- exits karyotype file
parser.add_argument('--chrid_list',help='chrid_list')
parser.add_argument('--karyotype',help='karyotype')
#---
parser.add_argument('--ideogram_conf',help='ideogram_conf',required=True)
parser.add_argument('--ticks_conf',help='ticks_conf',required=True)
parser.add_argument('--fusion_circos_conf',help='fusion_circos_conf',required=True)
parser.add_argument('--circos',help='circos',required=True)


argv = vars(parser.parse_args())
gmapx = argv['gmap']
fasta = argv['fasta']
gmap_idx = argv['gmap_idx']
gmap_idx_p = argv['gmap_idx_p']
gtf = argv['gtf']
STARx = argv['STAR']
clean_fq1 = argv['clean_fq1']
clean_fq2 = argv['clean_fq2']
genome = argv['genome']
abbr = argv['abbr']
num = argv['num']
chrid_list = argv['chrid_list']
karyotype = argv['karyotype']
ideogram_conf = argv['ideogram_conf']
ticks_conf = argv['ticks_conf']
fusion_circos_conf = argv['fusion_circos_conf']
circos = argv['circos']


indir = os.getcwd()
gmap_sum = indir+'/gmap_alignment.summary'
gmap_log = indir+'/gmap_alignment.log'
step1_result = indir+'/raw_result.xls'
step1_result_fa = indir+'/raw_result.fa'
if argv['karyotype'] != None:
	karyotype = 'karyotype = '+argv['karyotype']
if argv['chrid_list'] == None:
	chrid_list = indir+'/chrid.list'
step1_result_sam = indir+'/Aligned.out.sam'
fusion_result = indir+'/fusion_result.xls'
fusion_links = indir+'/fusion_links.txt'
circos_conf = indir+'/fusion_circos.conf'
#--------------------------

def gmap():
	cmd = '{gmapx} -D {gmap_idx} -d {gmap_idx_p} --cross-species --expand-offsets 1 -B 5 -K 8000 -f samse -t 4 -S {fasta} > {gmap_sum} 2> {gmap_log}'.format(gmapx=gmapx,gmap_idx=gmap_idx,gmap_idx_p=gmap_idx_p,fasta=fasta,gmap_sum=gmap_sum,gmap_log=gmap_log)
	print cmd
	os.system(cmd)

def get_align(align):
    dic = {}
    with open(align) as f:
        for line in f:
            if line[0] == '>':
                read_id = line[1:].strip()
                dic[read_id] = []
            m = re.search("  Path \d: query (\d+)\.\.(\d+).+genome (.+?):(.+)\.\.([\d,]+) ",line)
            l = re.search("Genomic pos:.+? \(([+-]) strand\)",line)
            n = re.search("Coverage: ([\d\.]+) ",line)
            if m:
                s = int(m.group(1))
                e = int(m.group(2))
                g = m.group(3)
                gs = int(m.group(4).replace(',',''))
                ge = int(m.group(5).replace(',',''))
                dic[read_id].append([s,e,g,gs,ge])
            if l:
                strand = l.group(1)
                dic[read_id][-1].append('')
                dic[read_id][-1].append(strand)
            if n:
                c = float(n.group(1))
                dic[read_id][-1][-2] = c
    return dic 

def fusion_filter(align_detail):
    if len(align_detail) == 1:
        return "no"
    total_c = 0
    chrs = []
    se = []
    gse = []
    dist = []
    for each in align_detail:
#each coverage > 10%
        if each[5] < 10:
            return "no"
        total_c += float(each[5])
        chrs.append(each[2])
        for gs_ge in gse:
#same chromosome distance > 100,000 bp
            dist.append(min(abs(gs_ge[0]-each[4]),abs(gs_ge[1]-each[3])))
        for s_e in se:
#unique mapping
            if len(set(range(each[0],each[1]))&set(range(s_e[0],s_e[1]))) > 1:
                return "no"
        gse.append([each[3],each[4]])
        se.append([each[0],each[1]])
#total coverage > 99%
    if total_c < 99:
        return "no"
#different chromosome
    if len(set(chrs)) > 1:
            return "yes"
    for each in dist:
        if each < 100000:
            return "no"
    return "yes"

def STAR():
    cmd1 = '{STAR} --runMode genomeGenerate --runThreadN 10 --genomeDir {indir} --genomeFastaFiles {step1_fasta}\n'.format(STAR=STARx,indir=indir,step1_fasta=step1_result_fa)
    cmd2= '{STAR} --runThreadN 10 --genomeDir {indir} --readFilesIn {clean_fq1} {clean_fq2} \n'.format(STAR=STARx,indir=indir,clean_fq1=clean_fq1,clean_fq2=clean_fq2)
    print cmd1
    print cmd2
    os.system(cmd1)
    os.system(cmd2)

def check_locus(id,juction_site,juction_site2):
    locus =  range(int(juction_site)-19,int(juction_site)+21)
    l = os.popen("grep '	%s	' %s"%(id,step1_result_sam))
    f = l.read().strip().split('\n')
    support_read = []
    for line in f[1:]:
        score = 0
        list = line.split('\t')
        leng = len(list[9])
        if int(list[3]) >(int(juction_site)+35-leng) and int(list[3])<(int(juction_site)-35) and int(list[3]) >(int(juction_site2)+35-leng) and int(list[3])<(int(juction_site2)-35):
            start = list[3]
            match = list[5]
            seq = list[9]
            tmp1 = re.split('[A-Z]',match)
            tmp2 = re.split('\d+',match)
            match_list = []
            for i in range(len(tmp1)-1):
                match_list.extend([tmp2[i+1]]*int(tmp1[i]))
            if match_list.count('M')> (leng-20) and (leng-match_list.count('M'))<30:
                for i in range(leng):
                    site = int(start)+i+1
                    if site in locus:
                        if match_list[i] == 'M':
                            score += 1
            if score >= 37:
                    support_read.append([id,list[0],start,match])
            if len(support_read) > 1:
                    return "yes"
    return "no"

def confirm():
        dic_chr = {}
        with open(chrid_list) as f:
                for line in f:
                        tmp = line.strip().split()
                        dic_chr[tmp[0]] = tmp[1]
	dic = {}
	dic_raw_result = {}
	f = open(step1_result,'r').readlines()
	for line in f[1:]:
		tmp = line.strip().split()
		id = tmp[0]
		dic_raw_result[id] = line
		if tmp[3] == tmp[10]:
			color = 'color=red\n'
		else:
			color  = 'color=green\n'
		if tmp[3] in dic_chr:
			chr_1 = dic_chr[tmp[3]]
		if tmp[10] in dic_chr:
			chr_2 = dic_chr[tmp[10]]
		dic[id] = [tmp[2],tmp[8],chr_1,tmp[5],tmp[5],chr_2,tmp[11],tmp[11],color]
	dic_fa = {}
	with open(step1_result_fa) as f:
		for line in f:
			if line[0] == '>':
				id = line[1:].strip()
				dic_fa[id] = ''
			else:
				dic_fa[id] += line.strip()
	
	result = 'transcriptID\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\n'
	circos = ''
	for id in dic.keys():
		validate = check_locus(id,dic[id][0],dic[id][1])
		if validate == "yes":
			result += dic_raw_result[id]
			circos += ('\t'.join(dic[id][2:])).replace("'","").replace(' ','')
	open(fusion_result,'w').write(result)
	open(fusion_links,'w').write(circos)

def novalidation():
    dic_chr = {}
    with open(chrid_list) as f:
        for line in f:
            tmp = line.strip().split()
            dic_chr[tmp[0]] = tmp[1]
    dic = {}
    dic_raw_result = {}
    f = open(step1_result,'r').readlines()
    for line in f[1:]:
        tmp = line.strip().split()
        id = tmp[0]
        dic_raw_result[id] = line
        if tmp[3] == tmp[10]:
            color = 'color=red\n'
        else:
            color  = 'color=green\n'
        if tmp[3] in dic_chr:
            chr_1 = dic_chr[tmp[3]]
        if tmp[10] in dic_chr:
            chr_2 = dic_chr[tmp[10]]
        dic[id] = [tmp[2],tmp[8],chr_1,tmp[5],tmp[5],chr_2,tmp[11],tmp[11],color]
    dic_fa = {}
    with open(step1_result_fa) as f:
        for line in f:
            if line[0] == '>':
                id = line[1:].strip()
                dic_fa[id] = ''
            else:
                dic_fa[id] += line.strip()

    result = 'transcriptID\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\n'
    circos = ''
    for id in dic.keys():
        result += dic_raw_result[id]
        circos += ('\t'.join(dic[id][2:])).replace("'","").replace(' ','')
    open(fusion_result,'w').write(result)
    open(fusion_links,'w').write(circos)

def plot():
	conf = open(fusion_circos_conf,'r').read()
	conf = conf.replace('karyotype',karyotype).replace('ideogram.conf',ideogram_conf).replace('ticks.conf',ticks_conf)
	open(circos_conf,'w').write(conf)
	cmd = 'perl {circos} -conf {circos_conf} '.format(circos=circos,circos_conf=circos_conf)
	print cmd
	os.system(cmd)

def generate_karyotype(fasta,abbr,num):
	global karyotype
	num = int(num)
	ends = []
	chrID_link = ''
	karyo = ''
	i = 0
	with open(fasta) as f:
		for line in f:
			if line[0]=='>':
				id = line[1:].strip()
				i += 1
				if i > num:
					break
				ends.append(0)
				chrID_link += id+'\t'+abbr+str(i)+'\n'
			else:
				ends[-1]+=len(line.strip())
	f = ''
	for i in range(1,num+1):
		colors = range(1,23)
		color = colors[(i-1)%22]
		karyo += 'chr - {abbr}{i} {i} 0 {end} chr{color}\n'.format(abbr=abbr,i=str(i),end=str(ends[i-1]),color=str(color))
	open('chrid.list','w').write(chrID_link)
	open('karyotype.txt','w').write(karyo)
	karyotype = 'karyotype = karyotype.txt'

def identify():
    dic_fa = {}
    with open(fasta) as f:
        for line in f:
            if line[0]==">":
                id = line[1:].strip()
                dic_fa[id] = line
            else:
                dic_fa[id] += line

    dic_align = get_align(gmap_sum)

    raw_result = 'transcriptID\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\tstart\tend\tchr\tchrstart\tchrend\tcoverage\tgene\n'
    raw_result_fa = ''
    geneModel = loadGeneModels(gtf,verbose=False)
    for id in dic_align.keys():
        filter = fusion_filter(dic_align[id])
        if filter == "yes":
            tmp = dic_align[id]
#            tmp[0].append(location2gene(tmp[0][2].strip().strip("'"),tmp[0][4],dic_location))
            genes = geneModel.getGenesInRange(tmp[0][2].strip().strip("'"),tmp[0][3],tmp[0][4],strand=tmp[0][-1])
            if len(genes) == 0:
                tmp[0][-1] = 'notsure'
            else:
                tmp[0][-1] = genes[0].id
#            tmp[1].append(location2gene(tmp[1][2].strip().strip("'"),tmp[1][3],dic_location))
            genes = geneModel.getGenesInRange(tmp[1][2].strip().strip("'"),tmp[1][3],tmp[1][4],strand=tmp[1][-1])
            if len(genes) == 0:
                tmp[-1][-1] = 'notsure'
            else:
                tmp[1][-1] = genes[0].id
            raw_result += id+'\t'
            raw_result += '\t'.join(['\t'.join(map(str,each)) for each in tmp])+'\n'
            raw_result_fa += dic_fa[id]
    dic_location = ''
    open(step1_result,'w').write(raw_result)
    open(step1_result_fa,'w').write(raw_result_fa)

'''
def gtf_location():
	mydic = {}
	l = os.popen('cut -f1 '+gtf)
	list = set(l.read().strip().split())
	for each in list:
		mydic[each] = {}
	with open(gtf) as f:
		for line in f:
			if '	gene	' in line:
				m = re.search('gene_id "(.+?)";',line)
				gene = m.group(1)
				tmp = line.split('\t')
				chr = tmp[0]
				s = int(tmp[3])
				e = int(tmp[4])
				for i in range(s,e,1000):
					mydic[chr][i]=gene
				mydic[chr][e]=gene
	return mydic
'''
'''
def location2gene(chr,num,dic):
	for i in range(2000):
		if num+i in dic[chr]:
			return dic[chr][num+i]
		if num-i in dic[chr]:
			return dic[chr][num-i]
	return 'notsure'
'''
#----------------------		

gmap()
identify()
if STARx != '' and clean_fq1 != '':
    STAR()
    if argv['karyotype'] == None:
	    generate_karyotype(genome,abbr,num)
    confirm()
else:
    if argv['karyotype'] == None:
        generate_karyotype(genome,abbr,num)
    novalidation()
plot()

