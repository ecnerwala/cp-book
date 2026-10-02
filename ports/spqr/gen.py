import random, sys
seed = int(sys.argv[1]); random.seed(seed)
NV = random.choice([1,2,3,4,5,6,8,10,12,20,40])
NE = random.randint(0, min(NV*NV, 60))
ternarize = random.randint(0,1)
edges = [(random.randrange(NV), random.randrange(NV)) for _ in range(NE)]
mode = seed % 3
vo = list(range(NV)); eo = list(range(NE))
if mode == 0: vo, eo = [], []
else:
    random.shuffle(vo); random.shuffle(eo)
    if mode == 2:
        vo = vo[:random.randint(0,NV)]; eo = eo[:random.randint(0,NE)]
print(NV, NE, ternarize)
for u,v in edges: print(u, v)
print(len(vo), *vo)
print(len(eo), *eo)
