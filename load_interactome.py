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
    counter_allnodes = 0
    counter_remainingnodes = 0
    for i in range(0,number_of_genes -1):
        node_id1 = int(disease_list[disease_id][i+4])
        for j in range(i+1,number_of_genes):
            counter_allnodes += 1
            node_id2 = int(disease_list[disease_id][j+4])

            if (node_id1 not in nx.nodes(G)) | (node_id2 not in nx.nodes(G)):
                continue

            if(nx.has_path(G,node_id1,node_id2)):
                counter_remainingnodes += 1
                average_shortest_path = average_shortest_path + nx.shortest_path_length(G,node_id1,node_id2)

    print ("number of nodes not skipped: ",counter_remainingnodes)
    print ("number of all nodes: ",counter_allnodes)
    average_shortest_path = average_shortest_path/float(counter_remainingnodes)
    return average_shortest_path




##########################################################################################################################
# code main body
#


# lets define our file path
f1 = open("DataS1_interactome_clean.tsv", "rb")
G  = lets_create_network_from_file(f1)
f1.close()


# generating a list of genes for each disease
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

nod  = len(disease_list)
avr  = np.zeros((nod,1))

for i in range(0,nod):
    avr[i] = compute_module_diameter(G,disease_list,i)
    print('#########################################')
    print("lets get the run down on disease:",disease_list[i][0], avr[i])
    print(' ')


########### Adding node attributes to each gene
# Reading Drug Bank drug targets; format: DrugBank ID, Name, Type, Uniprot ID, Uniprot Name
drugBank_file = np.genfromtxt("DrugBank_Uniprot_Data/drugbank_all_target_uniprot_links.txt",delimiter="\t",dtype=None)



#
