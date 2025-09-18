import sys
import gzip
#usage:
# python add_VAF_to_strelka2.py [strelka_vcf] [output_vcf] 
invcf=sys.argv[1]
outvcf=sys.argv[2]
tier2=sys.argv[3]

#https://github.com/Illumina/strelka/tree/master/docs/userGuide#somatic

if not (invcf.endswith(".gz")):
    print("WARNING! vcf files must be gzipped")
    exit()
elif "snvs" in invcf:

    indices=dict(AU=-4,CU=-3,GU=-2,TU=-1) #slicing index for XU

    formatcounter=0

    with gzip.open(invcf,'rt') as inv: #open input vcf
        with open(outvcf,'wt') as outv:       
            for line in inv:       
                if line.startswith('##') and (not line.startswith('##FORMAT')): #Info, contog vs is written directly
                    outv.write(line)
                elif line.startswith('##FORMAT') and formatcounter==0: #when we first encounter a format entry, I fist add my additional format info and then the rest of format entries are written
                    outv.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency of tier1 variants">\n')
                    outv.write('##FORMAT=<ID=AC,Number=1,Type=Integer,Description="alt allele count in tier 1">\n')
                    outv.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="ref,alt allele counts in tier 1">\n')
                    if tier2=="TRUE":
                        outv.write('##FORMAT=<ID=VAF2,Number=1,Type=Float,Description="Variant Allele Frequency of tier2 variants">\n')
                        outv.write('##FORMAT=<ID=AC2,Number=1,Type=Integer,Description="alt allele count in tier2 variants">\n')
                        outv.write('##FORMAT=<ID=AD2,Number=2,Type=Integer,Description="ref,alt allele count in tier2 variants">\n')
                    outv.write(line)
                    formatcounter+=1
                elif line.startswith('#'): #other header lines are unchanged
                    outv.write(line)
                else: #if they are not header, they are variants
                    line_list=line.strip().split("\t") # split the line into fields
                    ref=line_list[3] #ref allele
                    alt=line_list[4] # alt allele
                    tum=line_list[-1].split(':') # tumour info in the order of FMT 
                    normal=line_list[-2].split(':')# normal info in the order of FMT 
                    #### Tier 1 reads
                    # Retrieve allele counts A (alt) R (ref) for normal and tumour
                    AC1_normal=int(normal[indices[alt+'U']].split(",")[0]) #normal tier1 alternative allele count
                    AC1_tum=int(tum[indices[alt+'U']].split(",")[0])#tum tier1 alternative allele count
                    RC1_normal=int(normal[indices[ref+'U']].split(",")[0])#normal tier1 reference allele count
                    RC1_tum=int(tum[indices[ref+'U']].split(",")[0])#tum tier1 reference allele count
                    # Calculate VAF 
                    VAF1_normal= round(AC1_normal/(AC1_normal+RC1_normal),2) if AC1_normal+RC1_normal!=0 else 0 #tier1 VAF. if ref+alt is 0, VAF is 0
                    VAF1_tum=round(AC1_tum/(AC1_tum+RC1_tum),2) if AC1_tum+RC1_tum!=0 else 0
                    # Join AC in string for AD field
                    AD_tum=str(RC1_tum)+","+str(AC1_tum)
                    AD_normal=str(RC1_normal)+","+str(AC1_normal)
                    #### Tier 2 reads - optional
                    if tier2=="TRUE":
                        # Retrieve allele counts A (alt) R (ref) for normal and tumour
                        AC_all_normal=int(normal[indices[alt+'U']].split(",")[1])#normal tier2 alternative allele count
                        AC_all_tum=int(tum[indices[alt+'U']].split(",")[1])#tum tier2 alternative allele count
                        RC_all_normal=int(normal[indices[ref+'U']].split(",")[1])#normal tier2 reference allele count
                        RC_all_tum=int(tum[indices[ref+'U']].split(",")[1])#tum tier2 reference allele count
                        # Calculate VAF 
                        VAF2_normal=round(AC_all_normal/(AC_all_normal+RC_all_normal),2) if AC_all_normal+RC_all_normal!=0 else 0 #tier2 VAF
                        VAF2_tum=round(AC_all_tum/(AC_all_tum+RC_all_tum),2) if AC_all_tum+RC_all_tum!=0 else 0
                        # Join AC in string for AD field
                        AD2_tum=str(RC_all_tum)+","+str(AC_all_tum)
                        AD2_normal=str(RC_all_normal)+","+str(AC_all_normal)
                        # Create new lines 
                        fmt="DP:FDP:SDP:SUBDP:AU:CU:GU:TU:AC2:AD2:VAF2:AC:AD:VAF" #new FMT field 
                        new_tum=':'.join(tum)+':'+str(AC_all_tum)+":"+AD2_tum+":"+str(VAF2_tum)+":"+str(AC1_tum)+":"+AD_tum+":"+str(VAF1_tum) #new tumour field 
                        new_normal=':'.join(normal)+':'+str(AC_all_normal)+":"+AD2_normal+':'+str(VAF2_normal)+":"+str(AC1_normal)+":"+AD_normal+':'+str(VAF1_normal) # new normal field
                    # If tier2 not needed
                    else:
                        fmt="DP:FDP:SDP:SUBDP:AU:CU:GU:TU:AC:AD:VAF" #new FMT field 
                        new_tum=':'.join(tum)+':'+str(AC1_tum)+":"+AD_tum+":"+str(VAF1_tum) #new tumour field 
                        new_normal=':'.join(normal)+':'+str(AC1_normal)+":"+AD_normal+':'+str(VAF1_normal) # new normal field
                    newline="\t".join(line_list[0:8]+[fmt]+[new_normal]+[new_tum]) #new vcf line
                    outv.write(newline)
                    outv.write("\n")

elif "indels" in invcf:
    formatcounter=0
    with gzip.open(invcf,'rt') as inv: #open input vcf
        with open(outvcf,'wt') as outv:       
            for line in inv:       
                if line.startswith('##') and (not line.startswith('##FORMAT')): #Info, contog vs is written directly
                    outv.write(line)
                elif line.startswith('##FORMAT') and formatcounter==0: #when we first encounter a format entry, I fist add my additional format info and then the rest of format entries are written
                    outv.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency of tier1 variants">\n')
                    outv.write('##FORMAT=<ID=AC,Number=1,Type=Integer,Description="alt allele count in tier 1">\n')
                    outv.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="ref,alt allele counts in tier 1">\n')                    
                    if tier2=="TRUE":
                        outv.write('##FORMAT=<ID=VAF2,Number=1,Type=Float,Description="Variant Allele Frequency of tier2 variants">\n')
                        outv.write('##FORMAT=<ID=AC2,Number=1,Type=Integer,Description="alt allele count in tier2 variants">\n')
                        outv.write('##FORMAT=<ID=AD2,Number=2,Type=Integer,Description="ref,alt allele count in tier2 variants">\n')
                    outv.write(line)
                    formatcounter+=1
                elif line.startswith('#'): #other header lines are unchanged
                    outv.write(line)
                else: #if they are not header, they are variants
                    line_list=line.strip().split("\t") # split the line into fields
                    tum=line_list[-1].split(':') # tumour info in the order of FMT 
                    normal=line_list[-2].split(':')# normal info in the order of FMT 
                    #### Tier 1 reads
                    # Retrieve allele counts A (alt) R (ref) for normal and tumour
                    AC1_normal=int(normal[3].split(",")[0]) #normal tier1 alternative allele count
                    AC1_tum=int(tum[3].split(",")[0])#tum tier1 alternative allele count
                    RC1_normal=int(normal[2].split(",")[0])#normal tier1 reference allele count TAR
                    RC1_tum=int(tum[2].split(",")[0])#tum tier1 reference allele count
                    # Calculate VAF 
                    VAF1_normal= round(AC1_normal/(AC1_normal+RC1_normal),2) if AC1_normal+RC1_normal!=0 else 0 #tier1 VAF. if ref+alt is 0, VAF is 0
                    VAF1_tum=round(AC1_tum/(AC1_tum+RC1_tum),2) if AC1_tum+RC1_tum!=0 else 0
                    # Join AC in string for AD field
                    AD_tum=str(RC1_tum)+","+str(AC1_tum)
                    AD_normal=str(RC1_normal)+","+str(AC1_normal)
                    #### Tier 2 reads - optional
                    if tier2=="TRUE":
                        # Retrieve allele counts A (alt) R (ref) for normal and tumour
                        AC_all_normal=int(normal[3].split(",")[1])#normal tier2 alternative allele count
                        AC_all_tum=int(tum[3].split(",")[1])#tum tier2 alternative allele count TIR
                        RC_all_normal=int(normal[2].split(",")[1])#normal tier2 reference allele count
                        RC_all_tum=int(tum[2].split(",")[1])#tum tier2 reference allele count
                        # Calculate VAF 
                        VAF2_normal=round(AC_all_normal/(AC_all_normal+RC_all_normal),2) if AC_all_normal+RC_all_normal!=0 else 0 #tier2 VAF
                        VAF2_tum=round(AC_all_tum/(AC_all_tum+RC_all_tum),2) if AC_all_tum+RC_all_tum!=0 else 0
                        # Join AC in string for AD field
                        AD2_tum=str(RC_all_tum)+","+str(AC_all_tum)
                        AD2_normal=str(RC_all_normal)+","+str(AC_all_normal)
                        fmt="DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50:AC2:AD2:VAF2:AC:AD:VAF" #new FMT field 
                        new_tum=':'.join(tum)+':'+str(AC_all_tum)+":"+AD2_tum+":"+str(VAF2_tum)+":"+str(AC1_tum)+":"+AD_tum+":"+str(VAF1_tum) #new tumour field 
                        new_normal=':'.join(normal)+':'+str(AC_all_normal)+":"+AD2_normal+':'+str(VAF2_normal)+":"+str(AC1_normal)+":"+AD_normal+':'+str(VAF1_normal) # new normal field
                    else:
                        fmt="DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50:AC:AD:VAF" #new FMT field 
                        new_tum=':'.join(tum)+":"+str(AC1_tum)+":"+AD_tum+":"+str(VAF1_tum) #new tumour field 
                        new_normal=':'.join(normal)+':'+str(AC1_normal)+":"+AD_normal+':'+str(VAF1_normal) # new normal field
                    newline="\t".join(line_list[0:8]+[fmt]+[new_normal]+[new_tum]) #new vcf line
                    outv.write(newline)
                    outv.write("\n")
