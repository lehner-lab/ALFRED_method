import re,sys,os, glob
from string import *
import math, numpy, scipy, math
from numpy import array
from scipy import stats 
import pvalue_combine

##############
# This code is designed to detect LOH events of interested gene per one cancer patient using two different exome sequencing samples from normal and tumor.


#To run this code, below files need to prepare:

#1. ExAC variant file from http://exac.broadinstitute.org/downloads: ExAC.r0.3.sites.vep.vcf 
#2. gene2locus information file: UCSC-2014-10-30_hg19_refGene.txt
#3. cancer sequencing samples (please, see the example files from input directory), (1) normal sample from cancer patient, (2) tumor sample from cancer patient

#############


def germline2somatic_variant_mapping_LOHcalling (germline_sample, somatic_sample, chr_query, gene_query):
       
    #### step 1: collecting ExAC PASS variants
    ExAC_PASS = {}
    ExAC_PASS[chr_query] = []
    fexac = open('./input_folder/ExAC.r0.3.sites.vep.vcf','r') ## download from : http://exac.broadinstitute.org/downloads
    for line in fexac.xreadlines():
        line = line.strip()
        field = line.split('\t')
        if '##' in line:
            continue
        elif '#CHROM' in line:
            for i in range(len(field)):
                if field[i] == 'FILTER':
                    PASS_index = i
        else:
            chr_input = field[0]
            if chr_input!= chr_query:
                continue
            PASS_query = field[PASS_index]
            if PASS_query == 'PASS':
                variant_query = '%s-%s-%s-%s'%(field[0], field[1], field[3], field[4])
                ExAC_PASS[chr_input].append(variant_query)
    print len(ExAC_PASS[chr_query]), 'ExAC PASS', chr_query

    #### step 2: neighboring gene
    query_KB = 50*1000#to define neighboring variants. The definition of neighboring variants can be changed.
    flocus = open('./input_folder/UCSC-2014-10-30_hg19_refGene.txt','r')
    gene2locus = {}
    for line in flocus.xreadlines():
        line = line.strip()
        field = line.split('\t')
        gene_name = field[12]
        start_pos = int(field[4])
        end_pos = int(field[5])
        length = end_pos - start_pos
        chr_info = field[2]
        
        if gene_name not in gene2locus.keys():
            gene2locus[gene_name] = [[],[],[],[]]
            gene2locus[gene_name][0].append(length)
            gene2locus[gene_name][1].append(start_pos)
            gene2locus[gene_name][2].append(end_pos)
            gene2locus[gene_name][3].append(chr_info)
        else:
            gene2locus[gene_name][0].append(length)
            gene2locus[gene_name][1].append(start_pos)
            gene2locus[gene_name][2].append(end_pos)
    
    ### for finding neighboring genes
    gene2degree = {}
    degree2gene = {}
    gene_query_total = set([])
    
    for gids in gene_query:
        gene2degree[gids] = []
        gene_query_total.add(gids)
        
        length_query = numpy.median(gene2locus[gids][0])
        start_query = numpy.min(gene2locus[gids][1])
        end_query = numpy.max(gene2locus[gids][2])
        center_query = (start_query + end_query)/2.0
        chr_input = gene2locus[gids][3][0]
        possible_start = center_query - query_KB
        possible_end = center_query + query_KB

        for gene_cand in gene2locus.keys():
            chr_cand = gene2locus[gene_cand][3][0]
            if chr_input!= chr_cand:
                continue
            elif gene_cand == gids:
                continue
        
            length_cand = numpy.median(gene2locus[gene_cand][0])
            start_cand = numpy.min(gene2locus[gene_cand][1])
            end_cand = numpy.max(gene2locus[gene_cand][2])
        
            if end_cand > possible_start and end_cand < possible_end:
                gene2degree[gids].append(gene_cand)
                gene_query_total.add(gene_cand)
                if gene_cand in degree2gene.keys():
                    degree2gene[gene_cand].append(gids)
                else:
                    degree2gene[gene_cand] = [gids]
                    
            elif start_cand < possible_end and start_cand > possible_start:
                gene2degree[gids].append(gene_cand)
                gene_query_total.add(gene_cand)
                if gene_cand in degree2gene.keys():
                    degree2gene[gene_cand].append(gids)
                else:
                    degree2gene[gene_cand] = [gids]  
    
    print gene_query, 'gene_input_query'
    print gene_query_total, 'gene_neighboring_query' 
    
    ############################################################

    print germline_sample
    print somatic_sample
    
    fgermline = open('./input_folder/%s'%(germline_sample),'r')
    fsomatic = open('./input_folder/%s'%(somatic_sample), 'r')

    germline_variant = {}
    gene2variant = {}
    somatic_variant = {}
    ### collecting germline variants (all possible PASS variants)
    for line in fgermline.xreadlines():
        line = line.strip()
        field = line.split('\t')
        if '#' in line:
            if '#CHROM' in line:
                for i in range(len(field)):
                    if field[i] == 'FILTER':
                        filter_index = i
                    elif field[i] == 'INFO':
                        info_index = i
                    elif field[i] == 'Sample_index':
                        sample_index = i
                    elif field[i] == '#CHROM':
                        chr_index = i
        else:
            chr_input = field[chr_index]
            if chr_input!= chr_query:
                continue
            
            pos_query = field[chr_index+1]
            variant = '%s-%s-%s-%s'%(field[chr_index], field[chr_index+1], field[chr_index+3], field[chr_index+4])
            variant_index = '%s-%s-%s'%(field[chr_index], field[chr_index+1], field[chr_index+3])
            PASS_index = field[filter_index]
   
            gene_check = 0
            if PASS_index!= 'PASS':# to collect high-quality variant
                continue
            gene_annovar = field[info_index].split(';')#Annovar information
            
            for gindex in range(len(gene_annovar)):
                gene_info = gene_annovar[gindex]
                gids_info = gene_info.split('=')[0]
                if 'Gene.refGene' == gids_info:
                    gene_query_list = gene_info.split('=')[1].split(',')
                    for gids in gene_query_list:
                        if gids in gene_query_total:#collect variants of query genes and neighboring genes
                            gene_check = 1
                            
                elif 'Func.refGene' == gids_info:
                    exon = gene_info.split('=')[1]
                    
                elif 'ExonicFunc.refGene' == gids_info:
                    mut_type = gene_info.split('=')[1]
                    
                elif 'exac02' == gids_info:#ExAC MAF information
                    MAF = gene_info.split('=')[1]

            if gene_check == 0:
                continue                    
            elif variant not in ExAC_PASS[chr_input]:#to collect ExAC_PASS variants
                continue
            
            GQX_info = field[sample_index]
            genotype = GQX_info.split(':')[0]
            if genotype == '1/1':#to remove homozygous variant
                continue
            for gene_name in gene_query_list:
                gene2variant[gene_name] = [variant]
                germline_variant[variant] = ['%s\t%s\t%s\t%s\t%s'%(gene_name, mut_type, MAF, exon, GQX_info)]
                
    ### mapping germline variants from normal sample to matched somatic variants from tumor sample
    for line in fsomatic.xreadlines():
        line = line.strip()
        field = line.split('\t')
        if '#' in line:
            if '#CHROM' in line:
                for i in range(len(field)):
                    if field[i] == '#CHROM':
                        chr_index = i
                    elif field[i] == 'Sample_index':
                        sample_index = i
                    elif field[i] == 'FILTER':
                        filter_index = i
        else:
            chr_input = field[chr_index]
            if chr_input!= chr_query:
                continue
            profile = '%s-%s-%s-%s'%(field[chr_index], field[chr_index+1], field[chr_index+3], field[chr_index+4])            
            if profile in germline_variant.keys():
                somatic_GQX_info = field[sample_index]
                PASS_index = field[filter_index]
                somatic_variant[profile] = ['%s\t%s'%(PASS_index, somatic_GQX_info)]
             
    gene_info = {}
    for ids in gene_query:
        gene_info[ids] = [[]]

    for ids in germline_variant.keys():
        if ids in somatic_variant.keys():
            variant = ids
            ##############
            germline_mapping_info = germline_variant[ids][0].split('\t')

            gene = germline_mapping_info[0]
            mut_type = germline_mapping_info[1]
            MAF = germline_mapping_info[2]
            exon = germline_mapping_info[3]
            chr_input = ids.split('-')[0]

            ##############
            germline_GQX = germline_mapping_info[-1]
            germline_read_count = germline_GQX.split(':')[-1].split(',')
            germline_score1 = int(germline_read_count[0])
            germline_score2 = int(germline_read_count[1])

            somatic_mapping_info = somatic_variant[ids][0].split('\t')
            somatic_GQX = somatic_mapping_info[-1]            
            somatic_read_count = somatic_GQX.split(':')[-1].split(',')            
            somatic_score1 = int(somatic_read_count[0])
            somatic_score2 = int(somatic_read_count[1])

            if  float(somatic_score1 + somatic_score2) == 0.0:
                continue 
                    
            fisher_pvalue = stats.fisher_exact([[germline_score1, germline_score2], [somatic_score1, somatic_score2]])[1]
            variant_info = '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(variant,fisher_pvalue, mut_type, MAF, exon, germline_GQX, somatic_GQX)                 
            
            if gene in gene_info.keys():# collecting all possible variants of query gene
                gene_info[gene][0].append(variant_info)
            
            elif gene in degree2gene.keys():# collecting all possible variants of neighboring genes
                neighboring_list = degree2gene[gene]

                for friend_query in neighboring_list:
                    ###########################
                    start_pos_query = numpy.min(gene2locus[friend_query][1])
                    end_pos_query = numpy.max(gene2locus[friend_query][2])
                    
                    center_pos_query = (start_pos_query + end_pos_query)/2.0
                    possible_start = center_pos_query - query_KB
                    possible_end = center_pos_query + query_KB
                    ###########################
                    if position >= possible_start and position <= possible_end:
                        gene_info[friend_query][0].append(variant_info)
                            
                    elif position <= possible_end and position >= possible_start:
                        gene_info[friend_query][0].append(variant_info)
                        
    #################################
    effect_size_upper = 0.7 ## cut-off of effect size can be changed.
    effect_size_lower = 0.3
    gene_LOH = {}
    for gids in gene_info.keys():
        feature_list = gene_info[gids][0]
        gene_LOH[gids] = []
        
        #0: variant,1: fisher_pvalue, 2: mut_type, 3: MAF, 4: exon type, 5: germline_VAF, 6: somatic_VAF
        
        for sids in feature_list:
            mut = sids.split('\t')[2]
            exon = sids.split('\t')[4]

            germline_read = sids.split('\t')[-2]
            somatic_read = sids.split('\t')[-1]
            
            germline_info = germline_read.split(':')[-1].split(',')
            somatic_info = somatic_read.split(':')[-1].split(',')
            NormalVAF = int(germline_info[1]) / float(int(germline_info[0]) + int(germline_info[1]))
            TumorVAF = int(somatic_info[1]) / float(int(somatic_info[0]) + int(somatic_info[1]))
            fisher_pvalue = float(sids.split('\t')[1])
            if TumorVAF >= effect_size_upper or TumorVAF <= effect_size_lower:#effect size threshold to define LOH
                gene_LOH[gids].append(fisher_pvalue)
    
    fout = open('./output/LOH_mapping_output.txt', 'w')
    fout.write('Gene\tLOH_type\tNum_mapping_variant\n')

    for gene_name in gene_info.keys():
        pvalue_combination = pvalue_combine.combine_pvalues(gene_LOH[gene_name])[1]
        if pvalue_combination <= 0.05:
            fout.write('%s\t%s\t%s\n'%(gene_name, 'LOH', len(gene_LOH[gene_name])))
        else:
            fout.write('%s\t%s\t%s\n'%(gene_name, 'noLOH', len(gene_LOH[gene_name])))
    fout.close()
    print fout
    return 'LOH detection'
                    
#################################################

germline_sample = './input_folder/TCGA_GermlineSample_example.txt'
somatic_sample = './input_folder/TCGA_TumorSample_example.txt'
gene_query = ['BRCA1']
chr_query = '17'

print germline2somatic_variant_mapping (germline_sample, somatic_sample, chr_query, gene_query)

    
