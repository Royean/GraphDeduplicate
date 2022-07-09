
from sys import argv
import tqdm

in_graph = argv[1]
out_graph = argv[2]

L_adj = dict()
L_reorder = dict()
R_adj = dict()
R_reorder = dict()

with open(in_graph, 'r') as f:
    line = f.readline()
    while line:
        if line.startswith('%'):
            line = f.readline()
            continue
        info = line.split()
        u, v = info[0], info[1]
        if u not in L_reorder:
            # L_adj[u] = list()
            L_reorder[u] = len(L_reorder)
        if v not in R_reorder:
            R_reorder[v] = len(R_reorder)
        u,v = L_reorder[u], R_reorder[v]
        if u not in L_adj:
            L_adj[u] = set()
        if v not in R_adj:
            R_adj[v] = set()
        L_adj[u].add(v)
        R_adj[v].add(u)
        line = f.readline()
virtual_adj = R_adj if len(R_adj) > len(L_adj) else L_adj
with open(out_graph, 'w') as out:    
    sorted(virtual_adj)
    for k,v in tqdm.tqdm(virtual_adj.items()):
        out.write(str(k) + ' ')
        out.write(str(len(v)) + ' ')
        for u in v:
            out.write(str(u) + ' ')
        out.write(str(len(v)) + ' ')
        for u in v:
            out.write(str(u) + ' ')
        out.write("\n")