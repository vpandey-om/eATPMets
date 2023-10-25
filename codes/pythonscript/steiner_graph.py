import networkx as nx
from networkx.algorithms.approximation import steiner_tree
G = nx.Graph()

G.add_edge("a", "b", weight=1)
# G.add_edge("b", "c", weight=1)
G.add_edge("c", "f", weight=1)
G.add_edge("a", "d", weight=1)
# G.add_edge("d", "e", weight=1)
G.add_edge("e", "f", weight=1)

terminal_nodes=['c','f']
print(nx.is_connected(G))
print(nx.number_connected_components(G))
xx=[len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
import pdb; pdb.set_trace()
xx=steiner_tree(G, terminal_nodes)
print(xx.nodes)
# dist = lambda x, y: sum(abs(a - b) for a, b in zip(x, y))
# G1 = nx.random_internet_as_graph(9999)
# print('edges',len(G1.edges))

# terminal_nodes1=list(range(2,100,2))
# xx=steiner_tree(G1, terminal_nodes1)
# print(xx.nodes)
