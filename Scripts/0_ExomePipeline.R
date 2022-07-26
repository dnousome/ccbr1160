#############
#
#Run Exome pipeline on bam files (will need to unmap and remap to hg38)
#Run Updated pipeline? or CCBR pipeliner
##############


##Get the correct pairing 
readxl::read_xlsx("~/Desktop/WES_1160_parth.xlsx")

#git clone https://github.com/mtandon09/CCBR_GATK4_Exome_Seq_Pipeline.git

###Convert bams to fastq
#Tumor
tid=c("/data/SCLCgenomics/khanproc/atri_exo1/proc_7857901/processed_DATA/7857901/7857901/7857901_SB19_1389/7857901_SB19_1389.bwa.final.bam",
"/data/SCLCgenomics/khanproc/atri_exo1/proc_7944275/processed_DATA/7944275/7944275/7944275_SS19_1542/7944275_SS19_1542.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7975648/processed_DATA/7975648/7975648/sb_19_4740/sb_19_4740.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7997097/processed_DATA/7997097/7997097/sb_19_4812/sb_19_4812.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7996664/processed_DATA/7996664/7996664/sb_19_4900/sb_19_4900.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7984431/processed_DATA/7984431/7984431/sb_19_5365/sb_19_5365.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7955595/processed_DATA/7955595/7955595/sb_19_5931/sb_19_5931.bwa.final.bam")
out=basename(tid)
sapply(sprintf("cp %s ../bams/%s",tid,out),system)

#cd /data/CCBR/projects/ccbr1160/exome1/fastqs/normal
nid=c("/data/SCLCgenomics/khanproc/atri_exo1/proc_7857901/processed_DATA/7857901/7857901/7857901_N1D_E/7857901_N1D_E.bwa.final.bam",
"/data/SCLCgenomics/khanproc/atri_exo1/proc_7944275/processed_DATA/7944275/7944275/10_GUHA_MC_0438_001_N1D_E/10_GUHA_MC_0438_001_N1D_E.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7975648/processed_DATA/7975648/7975648/GUHA_MC_0522_1/GUHA_MC_0522_1.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7997097/processed_DATA/7997097/7997097/GUHA_MC_0600_1/GUHA_MC_0600_1.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7996664/processed_DATA/7996664/7996664/GUHA_MC_0573_1/GUHA_MC_0573_1.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7984431/processed_DATA/7984431/7984431/GUHA_MC_0566_1/GUHA_MC_0566_1.bwa.final.bam",
"/data/SCLCgenomics/khanproc/fred_exo1/proc_7955595/processed_DATA/7955595/7955595/GUHA_MC_0478_007_N1D_E/GUHA_MC_0478_007_N1D_E.bwa.final.bam")
out=basename(nid)
sapply(sprintf("cp %s ../normal/%s",nid,out),system)


##SAMTOOLS OR? GATK
##IF NO COLLAGE other takes up 22g
import os,glob,re
#bams=glob.glob("/data/desaipa/ccbr_analysis/bams/*.bam")
#bams=glob.glob("/vf/users/CCBR/projects/ccbr1160/exome1/fastqs/normal/*.bam")

b2fq=open("b2fq.swarm",'w')
for bam in bams:
  r=os.path.basename(bam)
  r1=re.sub(".bwa.final.bam",".R1.fastq.gz",r)
  r2=re.sub(".bwa.final.bam",".R2.fastq.gz",r)
  orphan=re.sub(".bwa.final.bam",".orphans.fastq.gz",r)
  c1=["gatk SamToFastq --INPUT",bam,"--FASTQ",r1,"--SECOND_END_FASTQ",r2,
          "--UNPAIRED_FASTQ",orphan,
          "--TMP_DIR","/data/CCBR/projects/ccbr1160/exome1/tmp",
          "-R","/data/CCBR_Pipeliner/Exome-seek/hg38/genome/Homo_sapiens_assembly38.fasta"]
  c1=' '.join(c1)
  b2fq.write(c1+"\n")

b2fq.close()

swarm -f b2fq.swarm -t 32 -g 48 --module GATK/4.2.4.1 --time=12:00:00


##RUN TUMOR ONLY mode
module purge
module load singularity snakemake

/data/CCBR/projects/ccbr1160/exome1/CCBR_GATK4_Exome_Seq_Pipeline/exome-seek run \
--dry-run \
--mode slurm \
--job-name wes_p0 \
--callers mutect2 vardict \
--input /data/CCBR/projects/ccbr1160/exome1/fastqs/*.R?.fastq.gz \
--output /data/CCBR/projects/ccbr1160/exome1/exomeseq \
--genome hg38 \
--targets /data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed 

/data/CCBR/projects/ccbr1160/exome1/CCBR_GATK4_Exome_Seq_Pipeline/exome-seek run \
--mode slurm \
--job-name wes_p0 \
--callers mutect2 vardict \
--input /data/CCBR/projects/ccbr1160/exome1/fastqs/*.R?.fastq.gz \
--output /data/CCBR/projects/ccbr1160/exome1/exomeseq \
--genome hg38 \
--targets /data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed 


##Run Germline TN mode as well!
tid=c("/data/SCLCgenomics/khanproc/atri_exo1/proc_7857901/processed_DATA/7857901/7857901/7857901_SB19_1389/7857901_SB19_1389.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/atri_exo1/proc_7944275/processed_DATA/7944275/7944275/7944275_SS19_1542/7944275_SS19_1542.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7975648/processed_DATA/7975648/7975648/sb_19_4740/sb_19_4740.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7997097/processed_DATA/7997097/7997097/sb_19_4812/sb_19_4812.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7996664/processed_DATA/7996664/7996664/sb_19_4900/sb_19_4900.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7984431/processed_DATA/7984431/7984431/sb_19_5365/sb_19_5365.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7955595/processed_DATA/7955595/7955595/sb_19_5931/sb_19_5931.bwa.final.bam")
nid=c("/data/SCLCgenomics/khanproc/atri_exo1/proc_7857901/processed_DATA/7857901/7857901/7857901_N1D_E/7857901_N1D_E.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/atri_exo1/proc_7944275/processed_DATA/7944275/7944275/10_GUHA_MC_0438_001_N1D_E/10_GUHA_MC_0438_001_N1D_E.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7975648/processed_DATA/7975648/7975648/GUHA_MC_0522_1/GUHA_MC_0522_1.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7997097/processed_DATA/7997097/7997097/GUHA_MC_0600_1/GUHA_MC_0600_1.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7996664/processed_DATA/7996664/7996664/GUHA_MC_0573_1/GUHA_MC_0573_1.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7984431/processed_DATA/7984431/7984431/GUHA_MC_0566_1/GUHA_MC_0566_1.bwa.final.bam",
      "/data/SCLCgenomics/khanproc/fred_exo1/proc_7955595/processed_DATA/7955595/7955595/GUHA_MC_0478_007_N1D_E/GUHA_MC_0478_007_N1D_E.bwa.final.bam")

n1=gsub(".bwa.final.bam","",basename(nid))
t1=gsub(".bwa.final.bam","",basename(tid))
p1=data.frame(Normal=n1,Tumor=t1)
readr::write_tsv(p1,"pairs.tsv")


/data/CCBR/projects/ccbr1160/exome1/CCBR_GATK4_Exome_Seq_Pipeline/exome-seek run \
--dry-run \
--mode slurm \
--cnv \
--job-name wes_p0_tn \
--pairs /vf/users/CCBR/projects/ccbr1160/exome1/fastqs/pairs.tsv \
--callers mutect2 strelka \
--input /data/CCBR/projects/ccbr1160/exome1/fastqs/fq/*.R?.fastq.gz \
--output /data/CCBR/projects/ccbr1160/exome1/TN_exomeseq \
--genome hg38 \
--targets /data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed 

/data/CCBR/projects/ccbr1160/exome1/CCBR_GATK4_Exome_Seq_Pipeline/exome-seek run \
--mode slurm \
--job-name wes_p0_tn \
--pairs /vf/users/CCBR/projects/ccbr1160/exome1/fastqs/pairs.tsv \
--callers mutect2 strelka \
--input /data/CCBR/projects/ccbr1160/exome1/fastqs/fq/*.R?.fastq.gz \
--output /data/CCBR/projects/ccbr1160/exome1/TN_exomeseq \
--genome hg38 \
--targets /data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed 



##Annotate the Germline Calls for the 7 WES samples
#python
import argparse,os,subprocess,glob,re,pandas,csv

input=os.getcwd()
mutin=glob.glob("*germline.vcf.gz")

data=pandas.read_csv('/vf/users/CCBR/projects/ccbr1160/exome1/fastqs/pairs.tsv',sep='\t')
norm=data["Normal"].values.tolist()

newout='anno.swarm'
open(newout, "w")
for mutect in mutin:
  fi=re.sub(".germline.vcf.gz","",mutect)
  if fi in norm:
    print fi
    vcf=re.sub(".vcf.gz",".vcf",mutect)
    outanno=fi+"_pass.maf"
    cmd=["table_annovar.pl",mutect,"$ANNOVAR_DATA/hg38 -thread 2 -buildver hg38 -vcfinput -nastring .",
     "-polish -out",fi,
     "-protocol gene,dbnsfp35a,dbnsfp41a,clinvar_20200419,cosmic92_coding,cosmic92_noncoding,exac03,exac03nontcga,gnomad30_genome -operation g,f,f,f,f,f,f,f,f"]
    cmd1=' '.join(cmd)
    with open(newout, "a") as outfile:
      outfile.write(cmd1+"\n")


