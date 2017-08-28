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
    print("There might be some missing values in the matrix, check output variable \'to map\'")
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
def compute_module_diameter(G,gene_list):

    number_of_genes = int(gene_list[1])

    list_of_nodes = []
    for i in range(0,number_of_genes):
        node_id1 = int(gene_list[i+4])

        if (node_id1  in nx.nodes(G)):
            list_of_nodes.append(node_id1)

    sub = nx.subgraph(G,list_of_nodes)

    connected  = max(nx.connected_component_subgraphs(sub), key = len)
    dist = connected.number_of_nodes()

    return dist


### this function takes in connected graph G, disease list, and disease id
### loops over all nodes in disease list[id] and computers shortest path between pairs of nodes
### returns module average shortest path between all node pairs
def compute_shortest_path_length(G,gene_list):
    number_of_genes = int(gene_list[1])
    average_shortest_path = 0.0
    genes_in_interactome = 0

    for i in range(0,number_of_genes):
        node_id1 = int(gene_list[i+4])


        if(node_id1 in nx.nodes(G)):
            min_dist = 1000000
            genes_in_interactome = genes_in_interactome + 1

            for j in range(0,number_of_genes):

                node_id2 = int(gene_list[j+4])
                if(node_id1 != node_id2):

                    if (node_id2 in nx.nodes(G)):

                        if(nx.has_path(G,node_id1,node_id2)):
                            dist = nx.shortest_path_length(G,node_id1,node_id2)

                            if(dist < min_dist):
                                min_dist = dist

            average_shortest_path = average_shortest_path + min_dist

    average_shortest_path = average_shortest_path/float(genes_in_interactome)
    return average_shortest_path


## This function will take geneID matrix and turn it into a dictionnary to be able to give nodes attributes
def dic_node_attributes(geneID_info_matrix, geneID_col, attribute_col):
    ## Removing geneIDs that do not existing in the network and those with MISSING entries
    geneID_info_matrix_trimmed = []
    for i in range(len(geneID_info_matrix)):
        if (geneID_info_matrix[i][geneID_col]=='MISSING'):
            next
        elif (int(geneID_info_matrix[i][geneID_col]) in nx.nodes(G)):
            geneID_info_matrix_trimmed.append(geneID_info_matrix[i])
    ## Creating a dictionnary of gene_IDs and attributes
    geneID_info_matrix_transpose = list(map(list, zip(*geneID_info_matrix_trimmed)))
    geneIDs_and_attribute_dic = dict(zip(list(map(int, geneID_info_matrix_transpose[geneID_col])), geneID_info_matrix_transpose[attribute_col]))
    return geneIDs_and_attribute_dic



##########################################################################################################################
# code main body
#


# lets define our file path
f1 = open("DataS1_interactome_clean_noMethod.tsv", "rb")
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

nod = len(disease_list)
avr  = np.zeros((nod,1))
dist = np.zeros((nod,1))

for i in range(0,len(disease_list)):
    if(disease_list[i][0] == "multiple sclerosis"):
        print('#########################################')
        print("lets get the run down on disease:",disease_list[i][0])

        avr[i] = compute_shortest_path_length(G,disease_list[i])
        dist[i] = compute_module_diameter(G,disease_list[i])
        print('shorest path length:',avr[i])
        print('module size:',dist[i])
        print(' ')

## Make the disease_list into an array of geneIDs and diseases
disease_gene_array = []
for i in range(0,len(disease_list)):
    disease_gene_array.extend([[x,disease_list[i][0]] for x in disease_list[i][4:]])

geneIDs_and_disease = dic_node_attributes(geneID_info_matrix=disease_gene_array, geneID_col=0, attribute_col=1)
nx.set_node_attributes(G, 'disease', geneIDs_and_disease)




########### Adding node attributes to each gene: protein names and Drug names
# Reading Drug Bank drug targets; format: DrugBank ID, Name, Type, Uniprot ID, Uniprot Name

## Mapping gene_IDs to DrugBank Uniprot IDs
drugBank_file = "DrugBank_Uniprot_Data/drugbank_all_target_uniprot_links.txt"
drugbankinfo = read_file(drugBank_file)


mapping_file = "DrugBank_Uniprot_Data/DrugBank_uniprot_Gene_Mapping.txt"
the_map = read_file(mapping_file)

drugbankinfo = mapping(to_map=drugbankinfo,the_map=the_map,col_tomap = 3,col_tomapto = 0)

# Writing mapped file to drugbank info to gene IDs
f2 = open("DrugBank_Uniprot_Data/drugbank_all_target_uniprot_GeneID.txt","w")
f2.write( '\t'.join(['DrugBankID','Name','Type','UniprotID', 'UniprotName','GeneID','\n']))
for i in range(len(drugbankinfo)):
    f2.write( '\t'.join([str(x) for x in drugbankinfo[i]]))
    f2.write('\n')

f2.close()

## Setting attributes to each node
geneIDs_and_drugs = dic_node_attributes(geneID_info_matrix=drugbankinfo, geneID_col=5, attribute_col=1)
geneIDs_and_proteinIDs = dic_node_attributes(geneID_info_matrix=drugbankinfo, geneID_col=5, attribute_col=4)
nx.set_node_attributes(G, 'drugs', geneIDs_and_drugs)
nx.set_node_attributes(G, 'proteinIDs', geneIDs_and_proteinIDs)

#
########### Adding node attributes to each gene: pathway name
########### Parsing XML files
#
import xml.etree.ElementTree as ET
tree = ET.parse('KEGG_Complement_Pathways/hsa04610.xml')
root = tree.getroot()
pathway_genes = []
for component in root:
    if (component.attrib['type']=='gene'):
        gene_nbres = [x.split("hsa:")[1] for x in component.attrib['name'].split(" ")]
        pathway_genes.extend(gene_nbres)

pathway_title = root.attrib['title']
print(pathway_title, "genes are: ", pathway_genes)

gene_pathway_map = [[x,pathway_title] for x in pathway_genes]

## Setting pathway attributes to each node
geneIDs_and_pathway = dic_node_attributes(geneID_info_matrix=gene_pathway_map, geneID_col=0, attribute_col=1)
nx.set_node_attributes(G, 'pathway', geneIDs_and_pathway)

'''
print('#########################################')
print("lets get the run down on the:",pathway_title)

pathway_list = [pathway_title,str(len(gene_pathway_map_transpose[0])),str(len(gene_pathway_map_transpose[0])),str(len(gene_pathway_map_transpose[0]))]
pathway_list.extend(gene_pathway_map_transpose[0])
avr_pathway = compute_shortest_path_length(G,pathway_list)
dist_pathway = compute_module_diameter(G,pathway_list)
print('shorest path length:',avr_pathway)
print('module size:',dist_pathway)
print(' ')
'''
node_relabel = dict(zip(G.nodes(),list(map(str, G.nodes()))))
H=nx.relabel_nodes(G,node_relabel)
nx.write_gml(H,"test.gml")
#
#
