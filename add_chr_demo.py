data = 'ASDF_1_CpG_filter.txt'
f = open(data,"r",encoding='utf-8')
datas = f.readlines()
datas = datas[1]
print(datas)
list = datas.split("\t")
print(list)
