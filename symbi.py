import json
import getopt
import sys
import copy

file = sys.argv[1]
print(file)
save = sys.argv[2]

left_dict = {}
# mid_dict = {}

with open(file) as f:
    dic = json.load(f)
    for l in dic["RECORDS"]:
        mid_node = l['movieid']
        left_node = l['actorid']
        if (left_dict.get(mid_node, None) == None):
            left_dict[mid_node] = list()
        left_dict[mid_node].append(left_node)

# print(left_dict)
tmp = dict()
count = 0

for key in left_dict.keys():
    l = left_dict[key]
    # print(left_dict[key])
    for idx,val in enumerate(l):
        tt = val
        if tmp.get(val, None) == None:
            tmp[val] = count
            count += 1
        l[idx] = tmp[val]

num = 0
right = dict()
left = dict()

for key in left_dict.keys():
    # l = list()
    left[num] = left_dict[key]

    for c in left_dict[key]:
        # print(c)
        if right.get(c,None) == None:
            # right[count]  = l
            right[c] = list()
        right[c].append(num)
        # print(right[c])
        # for l in left_dict[key]:
    # count += 1
    num += 1

# print(left)
# print(right)
with open(save,'w') as f:
    f.write(str(len(left)) + "\n")
    f.write(str(len(left)) + "\n")
    for key in left.keys():
        f.write(str(key) + " ")
        neighbor = set()
        # print(left[key])
        for r in left[key]:
            neighbor = neighbor | set(right[r])
            # print(neighbor)
        # left_w[key] = neighbor
        
        f.write(str(len(neighbor)) + " ")
        for x in neighbor:
            f.write(str(x) + " ")
        f.write("\n")