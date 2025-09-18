

#########################################
#### 1. Filter Strelka - PASS + AD 5 ####
#########################################


## Add VAF
cd variants/Strelka_raw
for vcf in *.vcf.gz
do
    out=$(echo $vcf | sed s/.vcf.gz/.VAF.vcf/g)
    python3.6 add_VAF_to_strelka2.py $vcf $out NOTIER2
done
    

### Filter 
mkdir Strelka_filtered
for vcf in *.VAF.vcf
do
    bcftools view -f PASS -i 'AC[1]>4' $vcf > Strelka_filtered/$vcf
done



###########################
#### 2. VEP annotation ####
###########################

curl -O https://ftp.ensembl.org/pub/release-89/variation/VEP/mus_musculus_vep_89_GRCm38.tar.gz
tar xzf mus_musculus_vep_89_GRCm38.tar.gz
mv mus_musculus mus_musculus_vep89
ln -s mus_musculus_vep89 mus_musculus

REF=ref/Mus_musculus.GRCm38.dna_sm.toplevel.fa

module purge
module load VEP/95.0-foss-2018b-Perl-5.28.0

# for VCF in *.VAF.vcf
while read VCF
do
    echo $VCF
    outvcf=$(echo $VCF | sed s/.VAF./.filtered.VEP./g)

    vep --input_file $VCF \
        --output_file ../VEP_annotated/$outvcf \
        --format vcf --vcf --tsl\
        --fasta $REF \
        --species mus_musculus \
        --pick --offline --cache --dir_cache dir/.vep --cache_version 89 \
        --plugin Frameshift --plugin Wildtype --dir_plugins dir/VEP_plugins
# done
done<vcfs_to_annotate.txt



### Concatenate snvs and indels 
cd ../VEP_annotated

module load BCFtools

for sample in $(ls *_snvs.filtered.VEP.vcf | cut -f1-5 -d_)
do
    bgzip ${sample}_somatic_snvs.filtered.VEP.vcf
    bgzip ${sample}_somatic_indels.filtered.VEP.vcf
    tabix ${sample}_somatic_snvs.filtered.VEP.vcf.gz
    tabix ${sample}_somatic_indels.filtered.VEP.vcf.gz
    bcftools concat -a -O z ${sample}_somatic_snvs.filtered.VEP.vcf.gz ${sample}_somatic_indels.filtered.VEP.vcf.gz > ${sample}.merged.annotated.vcf.gz
done


## Add GT field
for mergedVCF in ../VEP_annotated/*.merged.annotated.vcf.gz
do
    GTvcf=$(basename $mergedVCF | sed s/annotated.vcf.gz/annotated.GT.vcf.gz/g)
    python3.6 add_GT_to_strelka2.py $mergedVCF $GTvcf 
done



##################################
#### 3. Neoantigen prediction ####
##################################

### Neoantigen prediction


pathIn=neoag_input/
pathOut=neoantigens/
qu=$pathOut/qu
out=$pathOut/out
mkdir -p $pathIn $qu $out

for ANNOTATED_VCF in ${pathIn}*.merged.annotated.GT.vcf.gz
do
    sample=$(basename $ANNOTATED_VCF | sed s/StrelkaBP_//g | sed s/_TAIL.merged.annotated.GT.vcf.gz//g)

    jobID=pvac.${sample}
    jobName=${qu}/${jobID}.sh

    echo "#!/bin/bash" >  ${jobName}
    echo "module load Singularity/3.11.3" >> ${jobName}

    echo "singularity run -B /home /pvactools/pvactools_latest.sif \
        pvacseq run $ANNOTATED_VCF \
        TUMOR \
        H-2-Kb,H-2-Db \
        NetMHCpan \
        $pathOut/${sample} \
        -e1 8,9,10 \
        --iedb-install-directory /opt/iedb  \
        -b 500 -a sample_name --netmhc-stab " >> $jobName

    chmod +x ${jobName}
    submit.py -c ${jobName} -n ${jobID} -o ${out}/${jobID}.out -e ${out}/${jobID}.err -p cpu -t 24:00:00 -u 1 -m 4G -g ${qu}/${jobID}.sbatch
done



# Analyse and plot in R - neoantigen_analysis.R



################################
##### 4. SCNA burden check  ####
################################


# In R scna_analysis.R
