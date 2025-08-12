#!/opt/conda/envs/myenv/bin/python3.8
import pandas as pd
#import pandasql as pdsql
import glob
from itertools import chain
from pathlib import Path
import vcf2csv
import os
import shutil
import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys
import re
import yaml

def main(FILTER_DP,FILTER_VAF,SAMPLES,PATH,EFFECT_FILTER,bed,project,buildver,databases,exon_only,run_spliceai,run_annovar,THREADS,gene_list,email_cgi,token_cgi,cancertype,tumor_only, filters,oncoKB_token):#run_cgi):#,divide_output)
    run_cgi=False
    print(filters)
    if email_cgi and token_cgi:
        run_cgi=True
    elif (email_cgi and not token_cgi) or (not email_cgi and token_cgi):
        print("review if you entered all data to run CGI")
        return 1
    if buildver.lower()=="grch38":
        buildver="hg38"
    if buildver.lower()=="grch37":
        buildver="hg19"
    if buildver not in ["hg19","hg37","hg38"]:
        print("Please enter a valid genome version (hg19/hg37/hg38)")
        return 0
    if os.environ.get("SINGULARITY_NAME"):
        os.chdir("/home/")
    docker_scripts=os.path.realpath(os.path.dirname(__file__))+"/"
    docker_output="output_"+project ###############
    
    # Output files
    file_snpeff=docker_output+"/merged_variants_snpeff.txt"
    file_vep=docker_output+"/merged_variants_vep.txt"
    file_annovar=docker_output+'/ann.'+buildver+'_multianno.csv'
    file_annovar_exonic=docker_output+'/ann.'+buildver+'_exonic_multianno.csv'
    file_annovar_by_gene=docker_output+'/ann.'+buildver+'_by_gene_multianno.csv'
    file_cgi=docker_output+"/alterations.tsv"
    file_oncokb=docker_output+"/oncokb.csv"
    file_array= [file_snpeff,file_vep,file_annovar,file_annovar_exonic,file_annovar_by_gene,file_cgi,file_oncokb]
    if bed==None:
        bed=0
    #else:
    #    bed=docker_home+bed

    if gene_list:
    	if os.path.exists(gene_list):
        	gene=pd.read_csv(gene_list,header=None)
        	gene_list=gene[0].to_list()
    	else:
    		print ("Review if the gene list and its path are correct")
    		exit (0)

    if not docker_output.endswith("/"):
        docker_output+="/"
    if not os.path.isdir(docker_output):
        os.mkdir(docker_output)
    if not os.path.isfile(SAMPLES):
        print ("Review if the path for the sample file is correnct")
        exit(0)

    metasamples = pd.read_csv(SAMPLES, sep=",", header=0)
    metasamples.columns=metasamples.columns.str.upper()
    metasamples = metasamples[metasamples['LANE']==1]
    cod=[]
    for i,row in metasamples.iterrows():
        cod.append("_".join(row['FASTQ_1'].split("/")[-1].split("_")[0:2]))
    metasamples.insert(len(metasamples.columns),'Code',[i for i in cod])
    print('Loaded {} samples'.format(metasamples.shape[0]))
    print('Searching for variants in {}'.format(PATH))
    # FIRST parse VCF files and create a dictionary by variant key containing the variants info
    #########################
    # Extract data from VCF #
    #########################
    print("Extract vcf")

    VEP=[]
    patients_VEP=[]

    SNPEFF=[]
    patients_SNPEFF=[]
    # Collect vcf paths
    for index, row in metasamples.iterrows():
        sample_id = row['SAMPLE']
        sample_patient = row['PATIENT']
        #tumor_only = all(metasamples[metasamples['PATIENT'] == sample_patient]['STATUS'] == 1)
        if row['STATUS'] == 0:
            continue
        # Retrieve variants
        if tumor_only:
            for elem in Path(PATH).rglob('{}*filtered_snpEff.ann.vcf.gz'.format(sample_id)):
                SNPEFF.append(str(elem))
                patients_SNPEFF.append(sample_patient)
        else:
            
            gdna = metasamples[(metasamples['PATIENT'] == sample_patient)
                               & (metasamples['STATUS'] == 0)]['SAMPLE'].to_numpy()[0]
            for elem in chain(Path(PATH).rglob('{}_vs_{}*vcf.gz'.format(sample_id, gdna))):
                if "filtered" not in str(elem):
                    continue
                if "VEP" in str(elem).upper():
                    VEP.append(str(elem))
                    patients_vep.append(sample_patient)
                elif "SNPEFF" in str(elem).upper():
                    patients_SNPEFF.append(sample_patient)
                    SNPEFF.append(str(elem))
    if not os.path.exists(file_snpeff) and len(SNPEFF):
        vcf2csv.vcf2csv(SNPEFF,patients_SNPEFF,file_snpeff,FILTER_DP,FILTER_VAF,filters,EFFECT_FILTER,bed,gene_list,exon_only)
    if not os.path.exists(file_vep) and len(VEP):
        vcf2csv.vcf2csv(VEP,patients_VEP,file_vep,FILTER_DP,FILTER_VAF,filters,EFFECT_FILTER,bed,gene_list,exon_only)
    
    ###################
    # Execute Annovar 
    ###################
    if run_annovar and not os.path.isfile(file_annovar):
        print("RUN ANNOVAR")
        vcf2csv.imput_annovar(docker_output)
        vcf2csv.annovar(buildver,docker_output,docker_scripts,databases,THREADS,file_annovar)
        vcf2csv.clean_annovar(file_annovar,docker_output)
    #############################
    # Filter Exons from annovar #
    #############################
    if not os.path.isfile(file_annovar_exonic) and os.path.isfile(file_annovar) and exon_only:
        file_annovar=vcf2csv.annovar_filter_exons(file_annovar,file_annovar_exonic)

    ###############################
    # Filter by gene from annovar #
    ###############################
    if not os.path.isfile(file_annovar_by_gene) and gene_list and  os.path.isfile(file_annovar):
        file_annovar=vcf2csv.annovar_filter_by_gene(file_annovar,file_annovar_by_gene,gene_list)

    ####################
    # Execute CGI #
    ####################
    if run_cgi and not os.path.isfile(file_cgi):
        print ("Run CGI")
        vcf2csv.input_cgi(docker_output)
        vcf2csv.run_cgi(docker_output,email_cgi,token_cgi,buildver,cancertype)

    ####################
    # Execute SpliceAI #
    ####################
    spliceai_out=""
    spliceai_param=[]
    spliceai_OUTPUT=None
    if run_spliceai:
        print("Run Spliceai")    
        if  os.path.isfile(file_annovar): 
            vcf2csv.imput_spliceai(docker_output,pd.read_csv(file_annovar,usecols=[0,1,3,4]))
        else:
            vcf2csv.imput_spliceai(docker_output,pd.read_csv(file_snpeff,usecols=[1,2,3,4]))
        if buildver=="hg19":
            genome = 'grch37'
            url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz"
        else:
            genome = 'grch38'
            url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"

        downloaded="fasta/"+url.split("/")[-1]
        from_file=downloaded.replace(".gz","")
        fasta=genome+'.genome.fa' 
        fasta_path="fasta/"+fasta
        if not os.path.exists("fasta/"):
            os.mkdir("fasta")
        if not os.path.exists(fasta_path):
            cmd="wget "+url+" -P fasta/"
            vcf2csv.execute_cmd(cmd)
            cmd="gunzip "+downloaded
            vcf2csv.execute_cmd(cmd)
            cmd="mv "+from_file+" "+fasta_path
            vcf2csv.execute_cmd(cmd)


        for IMPUT in glob.glob(docker_output+"tmp_spliceai*"):
            spliceai_OUTPUT=IMPUT.replace("tmp_","")
            # Uncomend the following line to run in a cpu cluster
            vcf2csv.run_spliceai([IMPUT,spliceai_OUTPUT,genome,fasta])
            # comend the following two lines  to run in a cpu enviroment
            #spliceai_param.append([IMPUT,spliceai_OUTPUT,docker_home,genome,fasta])
        #######################################
        # Read SpliceAI output into DataFrame #
        #######################################
    
        spliceai_out=docker_output+"spliceai.csv"
        spliceai_df=pd.DataFrame()
        if not os.path.exists(spliceai_out):
            for spliceai in glob.glob(docker_output+"spliceai*"):
                if spliceai_df.empty:
                    spliceai_df = vcf2csv.parse_spliceai(spliceai)
                else:
                    spliceai_df=pd.concat([spliceai_df,vcf2csv.parse_spliceai(spliceai)])
            spliceai_df.to_csv(spliceai_out,sep="\t",index=False)
            del(spliceai_df)
   
        #if not 'chr' in spliceai_df['CHROM'][0]:
        #   spliceai_df['CHROM']='chr'+spliceai_df['CHROM']


    ##############
    # Run oncoKB #
    ##############

    if oncoKB_token:
        print ("Run OncoKb")
        vcf2csv.run_oncoKB(docker_output,oncoKB_token,buildver)

    
    ###############################################################
    # Read SnpEff and VEP results and marge with SplicaAI Results #
    ###############################################################
    print("Generate output")
    try:
        README=PATH+"README.txt"
        programs=PATH+"programs.txt"
        yml=PATH+"/pipeline_info/software_versions.yml"
        with open(yml, 'r') as file:
            prime_service = yaml.safe_load(file)
        tools=dict()
        for k in prime_service:
            for k2 in prime_service[k]:
                tools[k2]=prime_service[k][k2]
        #if run_spliceai:

        f=open(programs,"w")
        for k in tools:
            string=k+"\t"+tools[k]+"\n"
            f.write(string)
        f.close
    
        pipeline_info=PATH+"/pipeline_info/execution_report*html"
        html=glob.glob(pipeline_info)[0]
        f=open(html,"r")
        for line in f:
            if "nfcommand" in line and "class" in line:
                nextflow=line.split(">")[3].split("<")[0]+"\n"
        f.close()
        readme=open(README,"a")
        readme_txt=" ".join(sys.argv)+"\n"
        readme.write("nextflow comand:\n")
        readme.write(nextflow)
        readme.write("\nparse_variants.py command:\n")
        readme.write(readme_txt)
        readme.write("\nparse_variants parameters:\n")
        param=locals()
        for key in param:
            string=str(key)+" : "+str(param[key]).rstrip()+"\n"
            readme.write(string)
        readme.close()
    except:
        print ("I miss some files to generate the README")
    
    excel_OUTPUT=PATH+'merged_variants_'+project

    vcf2csv.write_excel(file_snpeff,file_vep,spliceai_out,file_annovar,file_cgi,file_oncokb,excel_OUTPUT)#,divide_output)

    print("output:",excel_OUTPUT+".xlsx")
    #shutil.rmtree(docker_output)

if __name__ == '__main__':

    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('--filter-dp', type=int, default=0, required=False,
                        help='Value for the DP filter (coverage). Default=0')

    parser.add_argument('--filter-vaf', type=float, default=0, required=False,
                        help='Value for the VAF filter in a percentage (Variant Allele Frequency). Default=0')

    parser.add_argument('--bed', type=str, default=0, required=False,
                        help='bed_file')

    parser.add_argument('--samples', type=str, default=None, required=False,
                        help='Path to the file with the samples info (metadata Sarek).')

    parser.add_argument('--path', type=str, default=None, required=False,
                        help='Path to the folder of the results from Sarek (Annotated folder).')


    parser.add_argument('--IncludeLowImpact', action='store_true', default=False, required=False,
                        help='Include Low impact variants Default=False')

    parser.add_argument('--gene_list', type=str, default=False, required=False,
                        help='list of genes of interest (one column text file)')

    parser.add_argument('--exon_only', action='store_true', default=False, required=False,
                        help='Filter out introns. Default=False')

    parser.add_argument('--run_spliceai', action='store_true', default=False, required=False,
                        help='the script runs spliceAi. Default=False.')
                        
    parser.add_argument('--run_annovar', action='store_true', default=False, required=False,
                        help='the script run Annovar. Default=False')
    parser.add_argument('--help_cancer_type_cgi', action='store_true', default=False, required=False,
                        help='List of cancer types to use the CGI annotator. Default=False')                    
    
    parser.add_argument('--annovar_databases', default=False, required=False,help='Indicate the databases to use by annovar to annotate.')


    parser.add_argument('--cancer_type', type=str, default="CANCER", required=False,
                        help='To run CGI.')                        
                    
    parser.add_argument('--email_cgi', type=str, default=False, required=False,
                        help='Mail that used to log into CGI. Default=False')
                        
    parser.add_argument('--token_cgi', type=str, default=False, required=False,
                        help='Token given by cgi to use the API. Default=False')    
                                        
    parser.add_argument('--genome', default=False, required=False,help='h19/hg38.')
    
    parser.add_argument('--thread', default=16, required=False,help='Number of treads. Default=16')


    parser.add_argument('--tumor_only', action='store_true', default=False, required=False,help='Set true for analize unpaired samples')
    parser.add_argument('--project', type=str, default=False, required=False,
                        help='Project name to genereate the output.')

    parser.add_argument('--filters', default="PASS", required=False,help='Indicate the tags you want to keep. THE FILTERS MUST BE WRITTEN SEPARATED BY COMMAS i.e: "PASS,clustered_events,haplotype"')
    parser.add_argument('--oncokb_token', default="", required=False,help='Add your oncoKB token')

                        
                        
    
    args = parser.parse_args()

    if args.help_cancer_type_cgi:

        types_of_cancer_cgi = {"CANCER":"Any cancer type", "HEMATO":"Hematologic malignancies", "SOLID":"Solid tumors", "LK":"Leukemia","LY":"Lymphoma","MAS":"Mastocytosis","MDS":"Myelodisplasic syndrome","MYMA":"Myeloma","CMP":"Chronic myeloprliferitive neoplasm","EOS":"Eosinophilia","TX":"Torax","US":"Urinary system","DS":"Digestive system","ES":"Endocrine system","FRS":"Female reproductory system","HNC":"Head an neck","MRS":"Male reproductory system","MSKV":"Musculoskeletal & vascular systems","NEU":"Neuroendocrine","NS":"Nervous system","SK":"Skin"}
        if args.help_cancer_type_cgi or not args.cancer_type in types_of_cancer_cgi.keys():
            print("\nPlease, type a valid command for cancer type. The options are the following:\n")

        for key, values in types_of_cancer_cgi.items():
           print(key,values,sep = "\t")
    
    elif not args.samples or not args.path or not args.genome or not args.project:
        print ("parse_variants.py: error: the following arguments are required: --samples, --path, --genome, --project")
    else:
        main(args.filter_dp,
         args.filter_vaf,
         args.samples,
         args.path,
         args.IncludeLowImpact,
         args.bed,
         args.project,
         args.genome,
         args.annovar_databases,
         args.exon_only,
         args.run_spliceai,
         args.run_annovar,
         args.thread,
         args.gene_list,
         args.email_cgi,
         args.token_cgi,
         args.cancer_type,
         args.tumor_only,
         args.filters,
         args.oncokb_token)#,args.run_cgi)#,args.divide_output)
