import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys


###
def lets_plot(graph):
    nx.draw(graph)
    plt.show()

def debug():
    sys.exit("debug mode")


##
##function to read TAB-SEPERATED files and if wanted convert to numpy array
def read_file(file,convert_to_nparray = False):
    f3 = open(file,encoding = 'utf-8',errors = 'replace')
    array = []
    for line in f3.readlines()[1:]:
        array_tmp = line.split('\t')
        array_tmp[-1] = array_tmp[-1].strip()
        array.append(array_tmp)
    #
    f3.close()
    if convert_to_nparray==True:
        import numpy as np
        array = np.asarray(array)
    return array

##

##This function is to map any new data to existing data; to_map is the data you wanna add new columns to
## the map is the new data that has a common column with to_map; col_tomap and col_tomapto is the common column between both matrices
## the to_map and the_map, respectively.
def mapping(to_map,the_map,col_tomap = 3,col_tomapto = 0):
    col_tomap = 3
    col_tomapto = 0
    for i in range(len(to_map)):
        ID_tomapto = to_map[i][col_tomap].strip()
        to_map[i].append('MISSING')
        #print(ID_tomapto)
        for j in range(len(the_map)):
            if (the_map[j][col_tomapto]==ID_tomapto):
                to_map[i][-1]=the_map[j][col_tomapto+1].strip()
    return to_map



### this functions takes in file path (list of edges, notes) and creates a networkX graph
### returns graph G
def lets_create_network_from_file(f1):
    G = nx.read_edgelist(f1, nodetype= int, data=(('method:',str),))
    nx.write_edgelist(G,'edgelist_attributes.edgelist')
    return G


### this function takes in connected graph G, disease list, and disease id
### loops over all nodes in disease list[id] and computers shortest path between pairs of nodes
### returns module average shortest path between all node pairs
### module size correspons to number of nodes in largest connected subgraph
def compute_module_diameter(G,disease_list,disease_id):

    number_of_genes = int(disease_list[disease_id][1])

    list_of_nodes = []
    for i in range(0,number_of_genes):
        node_id1 = int(disease_list[disease_id][i+4])

        if (node_id1  in nx.nodes(G)):
            list_of_nodes.append(node_id1)

    sub = nx.subgraph(G,list_of_nodes)

    connected  = max(nx.connected_component_subgraphs(sub), key = len)
    dist = connected.number_of_nodes()

    return dist


### this function takes in connected graph G, disease list, and disease id
### loops over all nodes in disease list[id] and computers shortest path between pairs of nodes
### returns module average shortest path between all node pairs
def compute_shortest_path_length(G,disease_list,disease_id):
    number_of_genes = int(disease_list[disease_id][1])
    average_shortest_path = 0.0
    genes_in_interactome = 0

    for i in range(0,number_of_genes):
        node_id1 = int(disease_list[disease_id][i+4])


        if(node_id1 in nx.nodes(G)):
            min_dist = 1000000
            genes_in_interactome = genes_in_interactome + 1

            for j in range(0,number_of_genes):

                node_id2 = int(disease_list[disease_id][j+4])
                if(node_id1 != node_id2):

                    if (node_id2 in nx.nodes(G)):

                        if(nx.has_path(G,node_id1,node_id2)):
                            dist = nx.shortest_path_length(G,node_id1,node_id2)

                            if(dist < min_dist):
                                min_dist = dist

            average_shortest_path = average_shortest_path + min_dist

    average_shortest_path = average_shortest_path/float(genes_in_interactome)
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
dist = np.zeros((nod,1))

for i in range(0,nod):
    if(disease_list[i][0] == "multiple sclerosis"):
        print('#########################################')
        print("lets get the run down on disease:",disease_list[i][0])

        avr[i] = compute_shortest_path_length(G,disease_list,i)
        dist[i] = compute_module_diameter(G,disease_list,i)
        print('shorest path length:',avr[i])
        print('module size:',dist[i])
        print(' ')



########### Adding node attributes to each gene
# Reading Drug Bank drug targets; format: DrugBank ID, Name, Type, Uniprot ID, Uniprot Name

## Mapping gene_IDs to DrugBank Uniprot IDs
drugBank_file = "DrugBank_Uniprot_Data/drugbank_all_target_uniprot_links.txt"
drugbankinfo = read_file(drugBank_file)


mapping_file = "DrugBank_Uniprot_Data/DrugBank_uniprot_Gene_Mapping.txt"
the_map = read_file(mapping_file)

drugbankinfo = mapping(to_map=drugbankinfo,the_map=the_map,col_tomap = 3,col_tomapto = 0)


f2 = open("DrugBank_Uniprot_Data/drugbank_all_target_uniprot_GeneID.txt","w")
f2.write( '\t'.join(['DrugBankID','Name','Type','UniprotID', 'UniprotName','GeneID','\n']))
for i in range(len(drugbankinfo)):
    f2.write( '\t'.join([str(x) for x in drugbankinfo[i]]))
    f2.write('\n')

f2.close()


#
