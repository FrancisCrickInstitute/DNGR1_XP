import sys
import gzip
#usage:
# python add_VAF_to_strelka2.py [strelka_vcf] [output_vcf]
inVCFfile=sys.argv[1]
outVCFfile=sys.argv[2]


if not (inVCFfile.endswith(".gz")):
    print("WARNING! vcf files must be gzipped")
    exit()
else:
    formatcounter=0
    with gzip.open(inVCFfile,'rt') as invcf: #open input vcf
        with gzip.open(outVCFfile,'wt') as outvcf:       
            for line in invcf:       
                if line.startswith('##') and (not line.startswith('##FORMAT')): #Info, contog vs is written directly
                    outvcf.write(line)
                elif line.startswith('##FORMAT') and formatcounter==0: #when we first encounter a format entry, I fist add my additional format info and then the rest of format entries are written
                    outvcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    outvcf.write(line)
                    formatcounter+=1
                elif line.startswith('#'): #other header lines are unchanged
                    outvcf.write(line)
                else: #if they are not header, they are variants
                    line_list=line.strip().split("\t") # split the line into fields
                    tum=line_list[-1] # tumour info in the order of FMT 
                    normal=line_list[-2]# normal info in the order of FMT 
                    fmt=line_list[-3]
                    new_fmt="GT:"+fmt
                    new_tum="0/1:"+tum #new tumour field 
                    new_normal="0/0:"+normal # new normal field
                    newline="\t".join(line_list[0:8]+[new_fmt]+[new_normal]+[new_tum])+"\n" #new vcf line
                    outvcf.write(newline)
