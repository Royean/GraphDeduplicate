import json

file = r'imdb_ijs_roles.json'
left_dict = {}
mid_dict = {}


with open(file) as f:
    dic = json.load(f)
    for l in dic["RECORDS"]:
        mid_node = l['l_partkey']
        left_node = l['l_orderkey']
        # if(mid_dict.get(left_node, None)) == None:
        #     mid_dict[left_node] = 1
        # else:
        #     print(left_node)
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



count = 0
# # print(sorted(intersec))
# # extract graph
with open('roles.tri','w') as f:
    for key in left_dict.keys():
        f.write(str(count) + " ")
        count += 1
        tmp_list = left_dict[key]
        length = len(tmp_list)
        f.write(str(length) + " ")
        for ele in tmp_list:
            f.write(str(ele) + " ")
        tmp_list = left_dict[key]
        f.write(str(len(tmp_list)) + " ")
        for ele in tmp_list:
            f.write(str(ele) + " ")
        f.write("\n")

# print(intersec)