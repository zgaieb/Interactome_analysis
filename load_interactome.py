import networkx as nx
import numpy as np



### this functions takes in file path (list of edges, notes) and creates a networkX graph
### returns graph G
def lets_create_network_from_file(f1):
    G = nx.read_edgelist(f1, nodetype= int, data=(('method:',str),))
    nx.write_edgelist(G,'edgelist_attributes.edgelist')
    return G



### this function takes in connected graph G, disease list, and disease id
### loops over all nodes in disease list[id] and computers shortest path between pairs of nodes
### returns module average shortest path between all node pairs
def compute_module_diameter(G,disease_list,disease_id):
    number_of_genes = int(disease_list[disease_id][1])
    average_shortest_path = 0.0
    counter = 0
    for i in range(0,number_of_genes -1):

        node_id1 = int(disease_list[disease_id][i+4])
        print(node_id1)
        for j in range(i+1,number_of_genes):
            node_id2 = int(disease_list[disease_id][j+4])

            print(i,j,node_id1,node_id2)
            counter += 1
            average_shortest_path = average_shortest_path + nx.shortest_path_length(G,node_id1,node_id2)


    average_shortest_path = average_shortest_path/float(counter)
    return average_shortest_path




##########################################################################################################################
# code main body
#


# lets define our file path
f1 = open("DataS1_interactome_clean.tsv", "rb")
G  = lets_create_network_from_file(f1)
f1.close()


# assign each node its associated diseases
disease_info_file = 'DataS2_disease_genes.tsv'
f2 = open(disease_info_file,"r")

disease_list = []
for line in f2.readlines()[1:]:
    elements = line.split('\t')
    disease_list_tmp = []
    disease_list_tmp.extend(elements[0:4])
    disease_list_tmp.extend(elements[4].strip().split(';'))
    if len(elements)==6:
        disease_list_tmp.extend(elements[5].strip().split(';'))
    #
    disease_list.append(disease_list_tmp)

print(disease_list[0])

nod  = len(disease_list)
avr  = np.zeros((nod,1))

for i in range(0,nod):
    avr[i] = compute_module_diameter(G,disease_list,i)
    print(avr[i])
