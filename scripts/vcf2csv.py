import shutil
import gzip
import sys
import pandas as pd
import numpy as np
import re
import itertools
import vcfpy
import os
import subprocess
#import multiprocessing as mp
import time
import requests
import zipfile
from random import seed
from random import randint
import json
import sys
import requests


def run_oncoKB(outdir,token,genome):
    INPUT=outdir+"tmp_anno.txt"
    if not os.path.exists(INPUT):
        imput_annovar(outdir)
    output=outdir+"oncokb.csv"
    if os.path.isfile(output): 
        return 1
    df = pd.read_csv(INPUT,sep = "\t", names=["CHROM","START","END","REF","ALT","ID"]) 
    oncgenicdf=pd.DataFrame(columns=["CHROM","POS","REF","ALT","OncoKB"])
    if genome=='hg19':
        genome='GRCh37'
    else:
        genome='GRCh38'
    header = {
        'Accept': 'application/json',
        'Authorization': token
    }
    server = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange?"
    for i,row in df.iterrows():
        ext='genomicLocation='+str(row["CHROM"])+'%2C'+str(row["START"])+'%2C'+str(row["END"])+"%2C"+row['REF']+"%2C"+row["ALT"]+'&referenceGenome='+genome+'&evidenceType=ONCOGENIC'
        r = requests.get(server+ext, headers=header)
        if not r.ok:
            #r.raise_for_status()
            oncgenicdf.loc[len(oncgenicdf.index)] = [row["CHROM"],row["START"],row['REF'],row["ALT"], "Could not be computed"]
            #sys.exit()
        else:
            oncoKg=r.json()
            oncgenicdf.loc[len(oncgenicdf.index)] = [row["CHROM"],row["START"],row['REF'],row["ALT"], oncoKg['oncogenic']] 
    oncgenicdf.to_csv(output,sep="\t",index=None)
    return 1

def run_spliceai(imput):#Function that makes spliceai works
    [IMPUT,spliceai_OUTPUT,genome,fasta]=imput
    cmd = 'spliceai -I {} -O {} -R fasta/{} -A {}'.format(IMPUT,spliceai_OUTPUT,fasta,genome)#Command line spliceAI, at the time spliceAI works by chromosome
    os.environ['PYTHON_EGG_CACHE'] = '/tmp/'
    if not os.path.exists(spliceai_OUTPUT):
        execute_cmd(cmd)
    os.unlink(IMPUT)

def parse_spliceai(filename):#Function that checks lines of the output of spliceai
    #reader = vcfpy.Reader.from_path(filename)
    header=["CHROM",
            "POS",
            "REF",
            "ALT",
            "Gene_Name",
            "DS_AG",
            "DS_AL",
            "DS_DG",
            "DS_DL",
            "DP_AG",
            "DP_AL",
            "DP_DG",
            "DP_DL"]
    row=[]
    vcf_in=open(filename,"r")
    get=0
    for record in vcf_in:#Iterates through the output of spliceAI
        record=record.rstrip()
        tmp=record.split("\t")
        if "CHROM" in record:#Avoid headers
            get = 1
            continue
        if get:#set each column to a variable
            chrom = tmp[0]
            start = tmp[1]
            ref = tmp[3]
            alt = tmp[4]
            [ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL] =  ["","","","","","","","","",""]
            if 'SpliceAI' not in record:
                [ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL] =  ["","","","","","","","","",""]#Adjust the fields to the format needed
            else:
                for splice in tmp[7].replace("SpliceAI=","").split(","):
                    tmp=[]
                    for tmp2replace in splice.split("|"):
                        if tmp2replace==".":
                            tmp.append("")
                        else:
                            tmp.append(tmp2replace)

                    [ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL] = tmp#Adjust the fields to the format needed


            row.append([
                chrom,
                start,
                ref,
                ALLELE,
                SYMBOL,
                DS_AG,
                DS_AL,
                DS_DG,
                DS_DL,
                DP_AG,
                DP_AL,
                DP_DG,
                DP_DL])
    vcf_in.close()
    return pd.DataFrame(row,columns=header)#Return the file created


def parse_ann_tumor_vs_wt(File,patient,tumor_vaf,tumor_dp,filtering,bed,gene_list):
    get=0
    if bed:
        bed=read_bed(bed)
    
    Header_1=['Patient','Sample','CHROM','POS','REF','ALT','FILTER','Reference_Reads','Alternative_Reads','Total_Reads','VAF_Alternative','Homozigocity']
    Header_2=['Patient','Tumor','Control','CHROM','POS','REF','ALT','FILTER',
            'Tumor_Reference_Reads','Tumor_Alternative_Reads','Tumor_Total_Reads','Tumor_VAF_Alternative','Tumor_Homozigocity',
            'Control_Reference_Reads','Control_Alternative_Reads','Control_Total_Reads','Control_VAF_Alternative','Control_Homozigocity']
    tmp=File.split("/")
    tmp=tmp[-1].split(".")
    tumor,control=tmp[0].split("_vs_")
    rows=[]
    samples=[]
    with gzip.open(File,'r') as fin:
        for line in fin:
            line=line.decode("utf-8")
            line=line.rstrip()
            if "##INFO=<ID=ANN" in line:
                t=line.split(",")
                match = re.match(r"^.*\'(.*)'.*$",t[-1])
                Header_ANN=match.group(1).split(" | ")[1:]
                Annotation_Impact_index=match.group(1).split(" | ").index('Annotation_Impact')
                Gene_Name_index=match.group(1).split(" | ")[1:].index('Gene_Name')
                Annotation_index=match.group(1).split(" | ")[1:].index('Annotation')
            if "##INFO=<ID=CSQ" in line:
                t=line.split(",")
                match = re.match(r"^.*Format: (.*)",t[-1])
                Header_ANN=Header_1+match.group(1).replace('">','').split('|')[1:]
                Annotation_Impact_index=match.group(1).replace('">','').split('|').index('Annotation_Impact')
                Gene_Name_index=match.group(1).split(" | ")[1:].index('SYMBOL')
                Annotation_index=match.group(1).split(" | ")[1:].index('Consequence')
            if "#CHROM" in line:
                get=1
                samples=line.split("\t")[9:]
                tumor_index=9
                control_index=9
                if control == samples[0]:  # Changed "in" by "==" for more robustness (Alba 27/06/24)
                    tumor_index=10
                else:
                    control_index=10
                Header=Header_1+Header_ANN
                if len(samples)==2:
                    Header=Header_2+Header_ANN
                continue
            if get:
                data=line.split("\t")
                CHROM=data[0]
                POS=data[1]
                ID=data[2]
                REF=data[3]
                ALT=data[4].split(",")
                QUAL=data[5]
                FILTER=data[6]
                accepted_filters=filtering.split(",")#All the possible filters are given to the parser
                if set(FILTER.split(";"))-set(accepted_filters):#Intersection between the filters given and the filter assigned to the variant
                	continue

                variants=[REF]+ALT
                FORMAT=data[8].split(":")
                hformat=dict()
                hformat['control']=dict()
                hformat['tumor']=dict()
                FORMAT2_control=data[control_index].split(":")
                FORMAT2_tumor=data[tumor_index].split(":")
                DP_i= len(tuple(itertools.takewhile(lambda x: "DP" not in x,FORMAT)))
                hformat['control']['DP']=FORMAT2_control[DP_i]
                hformat['tumor']['DP']=FORMAT2_tumor[DP_i]
                GT_control=FORMAT2_control[0]
                GT_tumor=FORMAT2_tumor[0]
                hformat['control']['AD']=dict()
                hformat['tumor']['AD']=dict()
                if "AD" in FORMAT:
                    AD_i= len(tuple(itertools.takewhile(lambda x: "AD" not in x,FORMAT)))
                    AD_control=FORMAT2_control[AD_i].split(",")
                    AD_tumor=FORMAT2_tumor[AD_i].split(",")
                    for v in range(0,len(variants)):
                        try:
                            hformat['control']['AD'][variants[v]]=AD_control[v]
                        except:
                            hformat['control']['AD'][variants[v]]='N/A'
                        try:
                            hformat['tumor']['AD'][variants[v]]=AD_tumor[v]
                        except:
                            hformat['tumor']['AD'][variants[v]]='N/A'
                t=data[7].split(';')
                first_idx =0
                if "ANN" in line:
                    first_idx = len(tuple(itertools.takewhile(lambda x: "ANN" not in x,t)))
                else:
                    first_idx = len(tuple(itertools.takewhile(lambda x: "CSQ" not in x,t)))
                ANN=t[first_idx].split("=")[1]
                for ann in ANN.split(','):
                    ALT=ann.split("|")[0]
                    Annotation_Impact=ann.split("|")[Annotation_Impact_index]
                    if ALT in hformat['tumor']['AD']:
                        ADalt_tumor=hformat['tumor']['AD'][ALT]
                    else:
                        ADalt_tumor=0

                    if REF in hformat['tumor']['AD']:
                        ADref_tumor=hformat['tumor']['AD'][REF]
                    else:
                        ADref_tumor=0

                    DP_tumor=hformat['tumor']['DP']
                    VAF_tumor=str(np.around((int(ADalt_tumor)/ int(DP_tumor)), 3)) if float(DP_tumor) > 0.0 else 0.0
                    if ALT in hformat['control']['AD']:
                        ADalt_control=hformat['control']['AD'][ALT]
                    else:
                        ADalt_control=0
                    
                    if REF in hformat['control']['AD']:
                        ADref_control=hformat['control']['AD'][REF]
                    else:
                        ADref_control=0

                    DP_control=hformat['control']['DP']

                    VAF_control=str(np.around((int(ADalt_control)/ int(DP_control)), 3)) if float(DP_control) > 0.0 else 0.0
                    
                    ANN =ann.split("|")[1:]
                    
                    GeneName=ANN[Gene_Name_index]
                    
                    if (filter_gene_list(gene_list,GeneName) and filter_vaf(VAF_tumor,tumor_vaf) and filter_DP(DP_tumor,tumor_dp) and filter_using_bed(CHROM,POS,bed,GeneName)):
                        line_out=[patient,tumor,control,CHROM,POS,REF,ALT,FILTER,
                            ADref_tumor,ADalt_tumor,DP_tumor,VAF_tumor,GT_tumor,
                            ADref_control,ADalt_control,DP_control,VAF_control,GT_control]+ann.split("|")[1:]
                        rows.append(line_out)
    df=pd.DataFrame(rows,columns=Header)
    return df

def parse_ann(file,patient,FILTER_DP,FILTER_VAF,filtering,bed=False,WESFilter=False,gene_list=False):#Function that parse the output of snpeff and  vep
    get=0
    

    Header=['Patient','Sample','CHROM','POS','REF','ALT','Caller','FILTER','ADref','ADalt','DPtotal','VAF','Homozigocity']#Fields of the file
    rows=[] #Rows that will be added to the dataframe that will be the output, if the line matches all the conditions given, it will be append to this list.
    fin=0 #Variable where the file will be set
    caller=0 #Variable where the caller will be set, mutect2, strelka, haplocaller, etc.
    ann_jocker=[]
    WES_Filter=[]
    if bed:#check if bed is given and read it if necessary
        bed=read_bed(bed)

    if "gz" in file:#In case the file is compressed, open with the function gzip
        fin=gzip.open(file,'r')
    else:
        fin=open(file,'r')
    nr=0
    for line in fin:#Check if it is in binary mode 
        if type(line)!=str:
            line=line.decode("utf-8")
        line=line.rstrip()
        if "##source=" in line:#After ##source the caller is given, so here it is rescued to use it afterwards
            t=line.split("=")
            caller=t[1]
        elif "##INFO=<ID=ANN" in line:#This pattern is the snpeff pattern, so as we know which annotator was used
            t=line.split(",")
            match = re.match(r"^.*\'(.*)'.*$",t[-1])
            Header=Header+match.group(1).split(" | ")[1:]
            for i in match.group(1).replace('">','').split('|'):
                ann_jocker.append("")
            Annotation_Impact=match.group(1).split(" | ")[1:].index('Annotation_Impact')
            Gene_Name_index=match.group(1).split(" | ")[1:].index('Gene_Name')
            Annotation_index=match.group(1).split(" | ")[1:].index('Annotation')
            WES_Filter=['stop_gained', 'missense_variant', '5_prime_UTR_premature_start_codon_gain_variant', 'disruptive_inframe_deletion', 'splice_acceptor_variant', 'splice_donor_variant', 'intron_variant', 'frameshift_variant', 'splice_region_variant', 'conservative_inframe_insertion', 'disruptive_inframe_insertion', 'conservative_inframe_deletion', '5_prime_UTR_variant', '3_prime_UTR_variant', 'stop_lost', 'non_coding_transcript_exon_variant', 'exon_loss_variant', 'start_lost', 'stop_retained_variant']
        elif "##INFO=<ID=CSQ" in line:#This pattern is the vep pattern, so as we know which annotator was used
            t=line.split(",")
            match = re.match(r"^.*Format: (.*)",t[-1])
            Header=Header+match.group(1).replace('">','').split('|')[1:]
            for i in match.group(1).replace('">','').split('|'):
                ann_jocker.append("")
            Annotation_Impact=match.group(1).split(" | ")[1:].index('IMPACT')
            Gene_Name_index=match.group(1).split(" | ")[1:].index('SYMBOL')
            Annotation_index=match.group(1).split(" | ")[1:].index('Consequence')
            WES_Filter=['splice_donor_variant','stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification','inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant','splice_region_variant','splice_donor_5th_base_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant','incomplete_terminal_codon_variant','start_retained_variant','stop_retained_variant','synonymous_variant','coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant','non_coding_transcript_exon_variant','intron_variant','NMD_transcript_variant','non_coding_transcript_variant','feature_elongation','feature_truncation']
        elif "#CHROM" in line:#set the samples in the vcf to a list.
            samples=line.split("\t")[9:]
        elif not line.startswith("#") and line!="":
            data=line.split("\t")
            CHROM=data[0]
            POS=data[1]
            ID=data[2]
            REF=data[3]
            ALT=data[4].split(",")
            QUAL=data[5]
            FILTER=data[6].split(";")
            accepted_filters=filtering.split(",")#All the possible filters are given to the parser
            if set(FILTER)-set(accepted_filters):#Intersection between the filters given and the filter assigned to the variant
                continue
            FILTER=";".join(FILTER)
            variants=[REF]+ALT
            FORMAT=data[8].split(":")
            i_format=9
            hformat=dict()#Dictionary of dictionaries (samples)
            for s in samples:
                hformat[s]=dict()
                FORMAT2=data[i_format].split(":")
                i_format+=1
                DP_i= len(tuple(itertools.takewhile(lambda x: "DP" not in x,FORMAT)))
                hformat[s]['DP']=FORMAT2[DP_i]#Set DP for the sample
                hformat[s]['GT']=FORMAT2[0]#Set GT for the sample --> Index GT = 0
                hformat[s]['AD']=dict()#Set AD for the sample
                AD_i= len(tuple(itertools.takewhile(lambda x: "AD" not in x,FORMAT)))
                try:
                    AD=FORMAT2[AD_i].split(",")#Set AD index
                except:
                    AD=np.nan
                for v in range(0,len(variants)):
                    try:
                        hformat[s]['AD'][variants[v]]=AD[v]#Set AD for the sample in the dictionary of the sample
                    except:
                        hformat[s]['AD'][variants[v]]='N/A'#Set AD for the sample in the dictionary of the sample
                t=data[7].split(';')
                first_idx =0
                ANN=[]
                if "ANN" in line:
                    first_idx = len(tuple(itertools.takewhile(lambda x: "ANN" not in x,t)))
                    ANN=t[first_idx].split("=")[1]#Set annotation line done by snpeff
                elif "CSQ" in line:
                    first_idx = len(tuple(itertools.takewhile(lambda x: "CSQ" not in x,t)))
                    ANN=t[first_idx].split("=")[1]#Set annotation line done by vep
                else:
                    ANN="|".join(ann_jocker)
                for ann in ANN.split(','):#Iterate through each of the different consequences
                    ALT=ann.split("|")[0]#Alternative allele for each cons
                    for s in samples:
                        DP=hformat[s]['DP']
                        if "strelka.genome" in file:#Strelka does not give AD and VAF
                            ADalt=np.nan
                            ADref=np.nan
                            VAF=np.nan
                        else:#Set AD and VAF for not strelka callers
                            ADalt=hformat[s]['AD'][ALT]
                            ADref=hformat[s]['AD'][REF]
                            VAF=np.around((int(ADalt)/ int(DP)) * 100, 3) if float(DP) > 0.0 else 0.0
                        GT=hformat[s]['GT']
                        ANN =ann.split("|")[1:]#List of the consequences/characteristics
                        GeneName=ANN[Gene_Name_index]
                        if ("strelka.genome" in file) or (filter_gene_list(gene_list,GeneName) and filter_vaf(VAF,FILTER_VAF) and filter_DP(DP,FILTER_DP) and filter_using_bed(CHROM,POS,bed,GeneName)):#Check the conditions given in the syntax
                            if not WESFilter or set(ANN[Annotation_index].split("&")).intersection(set(WES_Filter)):
                                line_out=[patient,s,CHROM,POS,REF,ALT,caller,FILTER,ADref,ADalt,DP,VAF,GT]+ANN
                                rows.append(line_out)#Append the lines that pass the filters
    df=pd.DataFrame(rows,columns=Header)#Updated each sample.
    
    return df

def filter_vaf(VAF,FILTER_VAF):#Function to filter by VAF
    if float(VAF)>=FILTER_VAF:
        return True
    return False

def filter_DP(DP,FILTER_DP):#Function to filter by DP
    if int(DP)>= FILTER_DP:
        return True
    return False

def filter_gene_list(gene_list,GeneName):#Function to filter by the list of the genes. If the list exists but the gene is not in the list return false, otherwise true (no list true and gene in list true)
    if gene_list and len(gene_list) > 0 and GeneName not in gene_list:
        return False
    return True

def filter_using_bed(CHROM,POS,bed,GeneName):
    if isinstance(bed, pd.DataFrame):    
        if len(bed.columns)>3:
            nameGeneCol=bed.columns[3]
            genelist=list(set(list(bed[nameGeneCol])))
        filter_gene=filter_gene_list(genelist,GeneName)
        inbed=False
        for i,row in bed.iterrows():
            if (CHROM==row['CHROM']) & (int(POS)>=row['START']) & (int(POS)>=row['END']):
                    inbed = True
        if filter_gene and inbed:
            return True
        return False
    return True

def read_bed(bed):#Function to read the bed file
    bed=pd.read_csv(bed,header=None,sep="\t")
    header=['CHROM','START','END']
    if len(bed.columns)>3:
        for col in range(3,len(bed.columns)):
            value_0=bed[col][0]
            if "NM" in value_0 or "ENSG" in value_0 or "ENSMUSG" in value_0:
                header.append('Feature_ID')
            elif "ENST" in value_0 or "ENSMUST" in value_0:
                header.append('ENST')
            elif "EXON" in value_0:
                header.append('EXON')
            elif "NM" in value_0:
                header.append('Feature_ID')
            else:
                header.append('Gene_Name')
    bed.columns=header
    return bed

def annovar_filter_exons(file_annovar,file_annovar_exonic):
    f_in=open(file_annovar,"r")
    f_out=open(file_annovar_exonic,"w")
    n=0
    #out=','.join(df.columns.to_list())+"\n"
    #f_out.write(out)
    for line in f_in:

        if "exonic" in line or n==0:
            f_out.write(line)
            n+=1
    f_in.close()
    f_out.close()
    return file_annovar_exonic
    
def annovar_filter_by_gene(file_annovar,file_annovar_by_gene,gene_list):
    f_in=open(file_annovar,'r')
    f_out=open(file_annovar_by_gene,"w")
    nl=0
    index=0
    for line in f_in:
        elements=line.split(",")
        if nl==0:
            nl=1
            index=elements.index('Gene.refGene')
            f_out.write(line)
            continue
        genes=elements[index].replace("\"","").split(";")
        if set(genes).intersection(gene_list):
            f_out.write(line)
    f_out.close()
    return file_annovar_by_gene
   
def annovar(buildver,docker_output,docker_scripts,databases,THREADS,file_annovar):
    if buildver=="hg37":
        buildver="hg19"
    IMPUT=docker_output+'tmp_anno.txt'#Create a temporal file for annovar
    ANNOVAR_PATH = docker_scripts+'table_annovar.pl'#Run annovar

    OUTPUT = docker_output+'/ann'#Prefix of the output
    #THREADS = 16
    if databases==False:
        databases=["refGene","avsnp150","cosmic92_coding","cosmic92_non_coding","gnomad_exome","MT_ensGene","refGeneWithVer","clinvar_20131105","dbnsfp41c","gnomad_genome","revel","clinvar_20220320","ljb26_all","refGeneVersion"]
    else:
        databases=databases.split(",")
        if "refGene" not in databases:
            databases.insert(0,"refGene")#RefGene is the most important database
    if not os.path.exists("./humandb/"):
        os.mkdir("humandb")
    listofdb=os.listdir("humandb/")
    for file in databases:
        subs=buildver+"_"+file
        #if len([x for x in listofdb if re.search(subs, x)])==0:  (Silencio temporal)
        #    cmd=docker_scripts+"annotate_variation.pl -buildver "+buildver+" -downdb --separate -webfrom annovar "+db+" humandb/"
        #    execute_cmd(cmd)
    f=["f" for db in range(1,len(databases))]
    f.insert(0,"g")

    cmd = '{} {} {} -buildver {} -thread {} -out {} -protocol {} -operation {} -nastring "" --nopolish -remove --checkfile --csvout'.format(
        ANNOVAR_PATH,
        IMPUT,
        'humandb',
        buildver,
        THREADS,
        OUTPUT,
        ','.join(databases),
        ','.join(f))#Be careful with the versions, hg19 is the only one that works.
    execute_cmd(cmd)
def generate_legend(col_dict,database,header):
    for i in header:
        if i in ['chr','CHROM','pos','POS','REF','ref1','ALT','alt1','ref','alt','Chr','Start','End','Ref','Alt']:
            continue
        col_dict[i]=database
    return col_dict


def remove_duplicated_columns(df):
    if '#Chr'  in list(df.columns):
        df.drop('#Chr',inplace=True,axis=1)
    for string in list(df.columns):
        m = re.search(r'\.\d', string)
        if m is not None:
            df.drop(string,inplace=True,axis=1)
    return df
def write_excel (file_snpeff,file_vep,spliceai,file_annovar,file_cgi,file_oncoKB,f_out):#,divide_output):#Function to create and configurate the output
    file_out=f_out+".xlsx"
    writer = pd.ExcelWriter(file_out, engine='xlsxwriter')
    writer.book.use_zip64()
    merged_snpeff=0
    merged_vep=0
    merged_annovar=0
    df_ann=0
    df_spliceai=0
    df_snpeff=0
    df_vep=0
    CGI=0
    df_oncokb=0
    field_origin=dict()
    if os.path.exists(file_annovar):
        df_ann=pd.read_csv(file_annovar)
        df_ann=remove_duplicated_columns(df_ann)
        field_origin=generate_legend(field_origin,"Annovar",list(df_ann.columns))
        c=list(df_ann.columns)
        c[0]="CHROM"
        c[1]="POS"
        c[3]="REF"
        c[4]="ALT"
        #c[-1]="Sample"
        df_ann.columns=c
        df_ann.drop(columns=['End'],inplace=True)
        #df_ann.drop(columns=['Gene.refGene'],inplace=True)
        colorder=list(df_ann.columns)
        colorder=colorder[-1:]+colorder[:-1]
        df_ann=df_ann[colorder]
    if os.path.exists(spliceai):#Create the file if it does not exist(No overwriting)
        df_spliceai=pd.read_csv(spliceai,sep="\t")
        field_origin=generate_legend(field_origin,"SpliceAI",list(df_spliceai.columns))
    if os.path.exists(file_snpeff):#Create the file if it does not exist(No overwriting)
        df_snpeff=pd.read_csv(file_snpeff,sep=',')
        field_origin=generate_legend(field_origin,"Snpeff",list(df_snpeff.columns))
    if os.path.exists(file_vep):
        df_vep=pd.read_csv(file_vep,sep=',')
        field_origin=generate_legend(field_origin,"VEP",list(df_vep.columns))
    if os.path.exists(file_oncoKB):
        df_oncokb=pd.read_csv(file_oncoKB,sep="\t")
        field_origin=generate_legend(field_origin,"OncoKB",list(df_oncokb.columns))
    if os.path.exists(file_cgi):
        CGI=pd.read_csv(file_cgi,sep="\t" , low_memory=False)
        field_origin=generate_legend(field_origin,"CGI",list(CGI.columns))
    '''
    if os.path.exists("./alterations.tsv"):#Create the file if it does not exist(No overwriting)
        Dataframe_cgi=pd.read_csv("./alterations.tsv", sep = "\t")
    else:
        with zipfile.ZipFile('CGI.zip', 'r') as zip_ref:
        	zip_ref.extractall()
        	Dataframe_cgi=pd.read_csv("./alterations.tsv", sep = "\t")      
    '''
    if isinstance(df_snpeff,pd.DataFrame) and isinstance(df_ann,pd.DataFrame): #Merge snpeff and annovar information

        df_snpeff[['CHROM', 'POS', 'REF', 'ALT']]=df_snpeff[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        df_ann[['CHROM', 'POS', 'REF', 'ALT']]=df_ann[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        merged=pd.merge(df_snpeff,df_ann,on=['CHROM', 'POS', 'REF', 'ALT'],how='right')
        df_snpeff=0
    elif isinstance(df_snpeff,pd.DataFrame):
        df_snpeff[['CHROM', 'POS', 'REF', 'ALT']]=df_snpeff[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        merged=df_snpeff#If there is no annovar file, the file will be only the snpeff file.
    
    if isinstance(CGI, pd.DataFrame):   

        CGI = CGI.rename(columns = {'chr':'CHROM','pos':'POS','REF':'ref1','ALT':'alt1','ref':'REF','alt':'ALT','ID':'Sample'})
        
        CGI[['CHROM', 'POS', 'REF', 'ALT']]=CGI[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        merged=merged.merge(CGI[['CHROM', 'POS', 'REF', 'ALT','CGI-Oncogenic Summary','CGI-Oncogenic Prediction']], on=['CHROM', 'POS', 'REF', 'ALT'],how='left')
    
    if isinstance(df_spliceai,pd.DataFrame):#Merge the merged file, merge the spliceAI to annovar and snpeff or only snpeff file.

        df_spliceai[['CHROM', 'POS', 'REF', 'ALT']]=df_spliceai[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        merged=pd.merge(merged,df_spliceai,on=['CHROM', 'POS', 'REF', 'REF', 'ALT','Gene_Name'],how='left')
    
    if isinstance(df_oncokb,pd.DataFrame):
        df_oncokb[['CHROM', 'POS', 'REF', 'ALT']]=df_oncokb[['CHROM', 'POS', 'REF', 'ALT']].astype(str)
        merged=pd.merge(merged,df_oncokb,on=['CHROM', 'POS', 'REF', 'REF', 'ALT'],how='left')
    merged['POS']=merged['POS'].astype(int)
    writer = format_excel_sheet(writer,merged,"results")#Call the function to configurate the excel
    if isinstance(df_vep,pd.DataFrame):#Remove possible duplicates
        merged_vep=merged_vep.drop_duplicates()
        merged_vep = merged_vep.loc[:,~merged_vep.columns.duplicated()].copy()
        writer = format_excel_sheet(writer,merged_vep ,"VEP")#Call the function for a VEP file
    legend=pd.DataFrame(columns=["Field","Origin"])
    for key in list(merged.columns):
        if  key in ["Patient","Sample"]:
            continue
        elif  key in ["CHROM","POS","REF","ALT","Caller","FILTER","ADref","ADalt","DPtotal","VAF","Homozigocity"]:
                legend.loc[len(legend.index)] = [key,"Caller"]
        else:
            legend.loc[len(legend.index)] = [key,field_origin[key]] 
    writer = format_excel_sheet(writer,legend,"Legend")#Call the function for a VEP file
    writer.close()
def format_excel_sheet(writer,df,name):
    df.to_excel(writer, sheet_name=name,startrow=1, header=False, index=False)#Design the file
    column_settings = [{'header': column} for column in df.columns]#Design the columns
    (max_row, max_col) = df.shape#Design the size
    workbook = writer.book
    worksheet = writer.sheets[name]#Select how many sheets and their titles
    header_format = workbook.add_format({'bold': True,'fg_color': '#00CCFF','border': 1})#Design the format of the header
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
    for col_num, value in enumerate(df.columns.values):
        worksheet.write(0, col_num, value, header_format)
        column_len = df[value].astype(str).str.len().max()
        # Setting the length if the column header is larger
        # than the max column value length
    column_len = max(column_len, len(value)) + 3
    worksheet.set_column(0, max_col - 1, 12)
    return writer

def imput_spliceai(outdir,spliceai_df):#Function to determine the input of spliceAI (at the time is splited by chromosome)

    tp='''##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n'''

    spliceai_df.columns=["CHROM","POS","REF","ALT"]
    spliceai_df=spliceai_df.replace("-",".")
    spliceai_df['ID']=["." for i,row in spliceai_df.iterrows()]
    spliceai_df=spliceai_df.reindex(["CHROM","POS","ID","REF","ALT"], axis=1)
    spliceai_df.drop_duplicates(inplace=True)
    #clean=outdir+"spliceai_CLEAN.csv"
    #spliceai_df.to_csv(clean)
    curent_chr=""
    f_spliceai=0
    for i,row in spliceai_df.iterrows():#Split by chromosome, create files
        if curent_chr!= row['CHROM']:
            if f_spliceai:
                f_spliceai.close()
            curent_chr = row['CHROM']
            splice_ai_tmp=outdir+"tmp_spliceai_"+row['CHROM']+".vcf"
            mode=1
            if os.path.exists(splice_ai_tmp):
                mode=0
            
            f_spliceai=open(splice_ai_tmp,"a")
            if mode:
                f_spliceai.write(tp)#Write headers
        towrite="\t".join([row["CHROM"],str(row["POS"]),row["ID"],row["REF"],row["ALT"]])+"\n"
        f_spliceai.write(towrite)
    f_spliceai.close()
    return 0

def execute_cmd(cmd):#All the commands are executed with this function
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        print(output)
        print(error)
        sys.exit(-1)
    return 0

def vcf2csv(files,patients,output,FILTER_DP,FILTER_VAF,filtering,EFFECT_FILTER,bed,gene_list,exon_only):#Function that converts the vcf file to a comma separated file
        head=0
        res=[i for i in files if "vs" in i]
        if res:
            files=res
        for i in range(len(files)):
            df=0
            if not "vs" in files[i]:
                df=parse_ann(files[i],patients[i],FILTER_DP,FILTER_VAF,filtering,bed,False,gene_list) #Be careful with the number of arguments here!!!
            else:
                df=parse_ann_tumor_vs_wt(files[i],patients[i],FILTER_VAF,FILTER_DP,filtering,bed,gene_list) #FILTER_VAF,FILTER_DP for tumor_vaf and tumor_dp
            if re.search('snpeff', files[i] , re.IGNORECASE) and EFFECT_FILTER==False:
                df=df[(df['Annotation_Impact']!='LOW')]
            elif re.search('vep', files[i] , re.IGNORECASE) and EFFECT_FILTER==False:
                df=df[(df['IMPACT']!='LOW')]
            #if re.search('snpeff', vcf , re.IGNORECASE) and exon_only==True:
            #    df=df[df['Annotation'].apply(lambda x: True if re.search('exon_variant', x) else False)]
            if not os.path.exists(output):
                df.to_csv(output,header=True, index=False)
            else:
                df.to_csv(output,header=False, mode='a', index=False)
        return 0

def imput_annovar(outdir):
    file_vep=outdir+"merged_variants_vep.txt"
    file_snpeff=outdir+"merged_variants_snpeff.txt"

    df=pd.DataFrame()
    if os.path.exists(file_vep):
        df=pd.read_csv(file_vep,usecols=["CHROM","POS","REF","ALT"])
    else:
        df=pd.read_csv(file_snpeff,usecols=["CHROM","POS","REF","ALT"])

    # df=df[list(df.columns[1:])+[df.columns[0]]]
    df.drop_duplicates(inplace=True)

    annovar_tmp=outdir+"tmp_anno.txt"
    df.fillna("-",inplace=True)
    df=df.replace(" ","-")
    df=df.replace("","-")
    ANN_I=df.iloc[:,0:5]
    END=[]
    for i,row in ANN_I.iterrows():
        END.append(int(row['POS'])+len(row['REF'])-1)
    ANN_I['END']=END
    
    ANN_I=ANN_I.reindex(['CHROM','POS','END','REF','ALT'],axis=1)
    ANN_I.drop_duplicates(inplace=True)
    if not os.path.exists(annovar_tmp):
        ANN_I.to_csv(annovar_tmp,sep="\t",index=False,header=None)
    else:
        ANN_I.to_csv(annovar_tmp,sep="\t",index=False,header=None,mode='a')
    return 0
def input_cgi(outdir):
    file_vep=outdir+"merged_variants_vep.txt"
    file_snpeff=outdir+"merged_variants_snpeff.txt"
    cgi_tmp=outdir+"tmp_cgi.txt"
    df=pd.DataFrame()
    if os.path.exists(file_vep):
        df=pd.read_csv(file_vep,usecols=["CHROM","POS","REF","ALT"])
    else:
        df=pd.read_csv(file_snpeff,usecols=["CHROM","POS","REF","ALT"])
    df=df[list(df.columns[1:])+[df.columns[0]]]
    df.drop_duplicates(inplace=True)
    df.fillna("-",inplace=True)
    df=df.replace(" ","-")
    df=df.replace("","-")
    CGI_I=df.iloc[:,0:5]
    END=[]
    for i,row in CGI_I.iterrows():
        END.append(int(row['POS'])+len(row['REF'])-1)
    CGI_I['END']=END
    CGI_I=CGI_I.reindex(['CHROM','POS','REF','ALT','Sample'],axis=1)
    CGI_I = CGI_I.rename(columns = {'CHROM':'chr','POS':'pos','REF':'ref','ALT':'alt','Sample':'ID'})
    CGI_I.to_csv(cgi_tmp,sep="\t",index=False)
    return 0

def run_cgi(outdir,email_cgi,token_cgi,genome,cancertype):
    headers = {'Authorization': email_cgi + " " + token_cgi}
    payload = {'cancer_type': cancertype, 'title': 'Annotation_CGI', 'reference': genome}
    r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',headers=headers,files={'mutations': open(outdir+"/tmp_cgi.txt", 'rb'),},data=payload)
    r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
    jobid=r.json()[-1]
    print("PROCESSING BY CGI jobid: "+jobid)
    while True:
        payload={'action':'logs'}
        r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/' + jobid, headers=headers, params=payload)
        if r.json()["status"] == "Done":
            break
        elif r.json()["status"] == "Error":
            print("Error in the CGI processing step, check your files or the status of the w    ebsite, then run the programm again")
            return 1
        else:
            time.sleep(randint(10,25))
    print ("Download Results from CGI")
    payload={'action':'download'}
    r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/'+jobid, headers=headers, params=payload)
    zip_file=outdir+"/cgi.zip"
    with open(zip_file, 'wb') as fd:
        fd.write(r._content)
    shutil.unpack_archive(zip_file, outdir)
    r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/'+jobid, headers=headers)
    r.json()
    return 0

def clean_annovar(file_annovar,tmp_dir):
    if not tmp_dir.endswith("/"):
        tmp_dir+="/"
    tmp=tmp_dir+"ann_tmp.csv"
    shutil.move(file_annovar,tmp)
    f=open(tmp,"r")
    f_out=open(file_annovar,"w")
    for line in f:
        line=line.strip().split(",")
        line=",".join(["" if c == "." else c for c in line])+"\n"
        f_out.write(line)
    f_out.close()
    f.close()
    os.remove(tmp)
