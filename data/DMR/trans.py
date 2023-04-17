def run(dirs_path):
    f = open(dirs_path , 'r')
    name = dirs_path.split('/')[-1].split('.')[0]
    w = open(name+'.dss.input.txt','w')
    w.write('chr'+'\t'+'pos'+'\t'+'N'+'\t'+'X'+'\n')
    for line in f:
        d = line.strip().split('\t')
        col_3 = int(d[-2])+int(d[-1])
        w.write(d[0]+'\t'+d[1]+'\t'+str(col_3)+'\t'+d[-1]+'\n')
    f.close()
    w.close()

run("/home/huchenyuan/ASD_project/data/methylation_information/ASDF_1_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/ASDF_2_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/ASDM_1_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/ASDM_2_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/CTF_1_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/CTF_2_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/CTM_1_trimmed_bismark_bt2.bismark.cov")
run("/home/huchenyuan/ASD_project/data/methylation_information/CTM_2_trimmed_bismark_bt2.bismark.cov")
