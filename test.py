import networkx as nx
import matplotlib.pyplot as plt

G = nx.path_graph(4)

seed_nodes = [1,2]

sub = nx.subgraph(G,seed_nodes)

comp = list(nx.connected_components(sub))

print(len(comp))
