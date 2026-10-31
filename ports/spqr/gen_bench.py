import random, sys
kind, NV, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); random.seed(seed)
edges = []
if kind == "random":   # random multigraph, NE = 2*NV -> mostly one big rigid (non-planar) component
    edges = [(random.randrange(NV), random.randrange(NV)) for _ in range(2 * NV)]
elif kind == "grid":   # planar grid, NV = W*W, one big rigid planar component
    W = int(NV ** 0.5); NV = W * W
    for i in range(W):
        for j in range(W):
            if i + 1 < W: edges.append((i * W + j, (i + 1) * W + j))
            if j + 1 < W: edges.append((i * W + j, i * W + j + 1))
elif kind == "sp":     # random series-parallel graph: lots of S/P nodes, planar
    edges = [(0, 1)]; nxt = 2
    while len(edges) < 2 * NV and nxt < NV:
        i = random.randrange(len(edges)); u, v = edges[i]
        if random.random() < 0.5: edges.append((u, v))
        else: edges[i] = (u, nxt); edges.append((nxt, v)); nxt += 1
    NV = nxt
elif kind == "tree":   # random tree + a few extra edges: bridges + small blocks
    edges = [(random.randrange(i), i) for i in range(1, NV)] + [(random.randrange(NV), random.randrange(NV)) for _ in range(NV // 10)]
print(NV, len(edges)); [print(u, v) for u, v in edges]
