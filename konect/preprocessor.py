# read the graph from cmd line
import sys

g = sys.argv[1]
save = sys.argv[2]

L_list = dict()
R_list = dict()

for line in open(g, 'r'):
    if line[0] == '%':
        continue
    edge = line.split()
    u,v = int(edge[0]), int(edge[1])
    if L_list.get(u, None) is None:
        L_list[u] = list()
    L_list[u].append(v)

    if R_list.get(v, None) is None:
        R_list[v] = list()
    R_list[v].append(u)

L_key2id = dict()
c = 0
for k in L_list.keys():
    if L_key2id.get(k, None) is None:
        L_key2id[k] = c
        c += 1

R_key2id = dict()
c = 0
for k in R_list.keys():
    if R_key2id.get(k, None) is None:
        R_key2id[k] = c
        c += 1

with open(save,'w') as out:
    out.write(str(len(L_key2id)) + "\n")
    out.write(str(len(R_key2id)) + "\n")
    for k in L_list.keys():
        out.write(str(L_key2id[k]) + " ")
        out.write(str(len(L_list[k])) + " ")
        for v in L_list[k]:
            out.write(str(R_key2id[v]) + " ")
        out.write("\n")