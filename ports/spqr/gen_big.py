import random, sys
seed = int(sys.argv[1]); random.seed(seed)
NV = random.choice([100, 300, 1000])
NE = random.randint(NV // 2, 3 * NV)
ternarize = seed % 2
edges = [(random.randrange(NV), random.randrange(NV)) for _ in range(NE)]
vo = list(range(NV)); eo = list(range(NE)); random.shuffle(vo); random.shuffle(eo)
vo = vo[:random.randint(0, NV)]; eo = eo[:random.randint(0, NE)]
print(NV, NE, ternarize)
for u, v in edges: print(u, v)
print(len(vo), *vo)
print(len(eo), *eo)
