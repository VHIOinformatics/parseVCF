# parseVCF

parseVCF takes a vcf file resulting from Mutect2 variant calling and annotated with SnpEff, which contains one variant per line, and converts it into Excel format with one line per annotation. It also allows the filtration and further annotation of variants.

The vcf file might come from somatic paired variant calling or from tumor-only variant calling; this must be specified with the --tumor_only argument.

The program is able to run ANNOVAR, SpliceAI, CGI and OncoKB to add extra annotations.


## Usage

parse_variants.py [-h] [--filter-dp FILTER_DP] [--filter-vaf FILTER_VAF] [--bed BED] \
		  [--samples SAMPLES] [--path PATH] [--IncludeLowImpact] [--gene_list GENE_LIST] \
                  [--exon_only] [--run_spliceai] [--run_annovar] [--help_cancer_type_cgi] \
                  [--annovar_databases ANNOVAR_DATABASES] [--cancer_type CANCER_TYPE] [--email_cgi EMAIL_CGI] \
                  [--token_cgi TOKEN_CGI] [--genome GENOME] [--thread THREAD] [--tumor_only] [--project PROJECT] \
		  [--filters FILTERS] [--oncokb_token ONCOKB_TOKEN]

optional arguments:
  -h, --help              show this help message and exit \
  --filter-dp FILTER_DP   Value for the DP filter (coverage). Default=0 \
  --filter-vaf FILTER_VAF Value for the VAF filter in a percentage (Variant Allele Frequency). Default=0 \
  --bed BED               bed_file \
  --samples SAMPLES       Path to the file with the samples info (metadata Sarek). \
  --path PATH             Path to the folder of the results from Sarek (Annotated folder). \
  --IncludeLowImpact      Include Low impact variants Default=False \
  --gene_list GENE_LIST   list of genes of interest (one column text file) \
  --exon_only             Filter out introns. Default=False \
  --run_spliceai          the script runs spliceAi. Default=False. \
  --run_annovar           the script run Annovar. Default=False \
  --help_cancer_type_cgi  List of cancer types to use the CGI annotator. Default=False \
  --annovar_databases ANNOVAR_DATABASES Indicate the databases to use by annovar to annotate. \
  --cancer_type CANCER_TYPE To run CGI. \
  --email_cgi EMAIL_CGI Mail that used to log into CGI. Default=False \
  --token_cgi TOKEN_CGI Token given by cgi to use the API. Default=False \
  --genome GENOME       h19/hg38. \
  --thread THREAD       Number of treads. Default=16 \
  --tumor_only          Set true for analize unpaired samples \
  --project PROJECT     Project name to genereate the output. \
  --filters FILTERS     Indicate the tags you want to keep. THE FILTERS MUST BE WRITTEN SEPARATED BY COMMAS i.e: "PASS,clustered_events,haplotype" \
  --oncokb_token ONCOKB_TOKEN Add your oncoKB token \


## Singularity image

There is a Singularity image specifically created to run the program: vcf.sif


To check the help of the main script:

> singularity run -H $PWD:/home/ vcf.sif python3 parseVCF/scripts/parse_variants.py --help


Example to run the program on tumor-only results on genome hg38 filtering variants with VAF < 0.05, keeping only variants with FILTER "PASS" or "germline" and keeping all values of the SnpEff Annotation_Impact field:

> singularity run -H $PWD:/home/ -B $sarek_results_dir:/home/results_dir/ vcf.sif python3 parseVCF/scripts/parse_variants.py --samples $samples --path results_dir/ --genome hg38 --filter-vaf 5 --filters "PASS,germline" --IncludeLowImpact --project $SLURM_JOB_NAME --tumor_only



## For developers

If you find a bug or there's a new feature that you would like to incorporate, open a GitHub issue so everyone is aware of it. It can be used as a communication channel to let others know that you are working on something.

To make modifications follow these steps:

1- On your local copy of the repository create a branch with your name and switch to it.

> git checkout -b <your_branch>


2- Make modifications and commit changes as usual.

> git add .

> git commit .

> git push


3- Test the modifications by running the program on your branch using the test data.


4- If everything works as expected, merge your branch to the development branch. First pull potential changes incorporated in the remote repository by other developers.

> git checkout development

> git pull

> git merge <your_branch>

> git push


5- Test the development branch, which might contain changes from multiple developers.


6- If everything works as expected, update program help and documentation and merge the development branch to the main branch.

> git checkout main

> git merge development

> git push


7- In the main branch create an annotated tag with the name of the program version and a description of the changes (it works like the commit command; if you don't specify a tagging message with -m, Git launches an editor to type the message). Push it to the remote repository.

> git tag -a v1.0.0

> git push origin --tags


8- To keep a more detailed description of the new version you can also create a "release" based on the tag you just added. This has to be done in GitHub (see https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)


