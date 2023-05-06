import sys
from Bio import SeqIO
import gzip

fq1 = "../ASD/ASDF_1.fastq"
fq2 = "../trim_data/ASDF_1_trimmed.fq"

for f in [fq1, fq2] :
    reads_count = 0
    gc_base = 0
    n_base = 0
    all_base = 0
    all_qual = 0
    q20_qual = 0
    q30_qual = 0
    length = []
    record = SeqIO.parse(f,'fastq')
    for lane in record :
        ## 序列数
        reads_count+=1
        ### 碱基统计
        gc_base+=len([i for i in lane.seq if i=='G' or i=='C'])
        n_base+=len([i for i in lane.seq if i=='N'])
        all_base+=(len(lane.seq))
        length.append(len(lane.seq))
        ### 质量值统计
        qual = lane.letter_annotations['phred_quality']
        all_qual+=(len(qual))
        q20_qual+=len([q for q in qual if q >=20])
        q30_qual+=len([q for q in qual if q >=30])
    #计算比例
    q20 = round(q20_qual/all_qual,4) * 100
    q30 = round(q30_qual/all_qual,4) * 100
    gc = round(gc_base/all_base,4) * 100
    n = round(n_base/all_base,4) * 1000000
    print(reads_count,all_base,gc,q20,q30,n,min(length),max(length))
