import json
import getopt
import sys

file = sys.argv[1]
print(file)
save = sys.argv[2]

left_dict = {}
mid_dict = {}

with open(file) as f:
    dic = json.load(f)
    for l in dic["RECORDS"]:
        mid_node = l['cited_paper_id']
        left_node = l['citing_paper_id']
        if (left_dict.get(mid_node, None) == None):
            left_dict[mid_node] = list()
        left_dict[mid_node].append(left_node)

tmp = dict()
count = 0

for key in left_dict.keys():
    l = left_dict[key]
    for idx,val in enumerate(l):
        tt = val
        if tmp.get(val, None) == None:
            tmp[val] = count
            count += 1
        l[idx] = tmp[val]

count1 = 0

with open(save,'w') as f:
    f.write(str(len(left_dict)) + "\n")
    f.write(str(count) + "\n")
    S = set(range(0, count))
    for key in left_dict.keys():
        f.write(str(count1) + " ")
        count1 += 1
        tmp_list = left_dict[key]
        length = count - len(tmp_list)
        f.write(str(length) + " ")
        for ele in S.difference(set(tmp_list)):
            f.write(str(ele) + " ")
        f.write("\n")
