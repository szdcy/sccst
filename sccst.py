import csv
import numpy as np
import re


def read_total_gene_in_human():
    # Create an empty list to store the gene information
    total_gene = []
    # total_gene = [[gene 1 name, gene 1 description, [gene 1 nickname 1, ...]], ...]

    # Open the CSV file containing the gene information
    with open("Homo_sapiens_gene_info.csv", newline='') as f:
        # Read the CSV file using the csv module
        reader = csv.reader(f)
        # Convert the reader object to a list
        raw_data = list(reader)
        # Iterate over each gene in the raw data
        for gene in raw_data:
            # Check if the gene has aliases
            if gene[2] != "-":
                # If the gene has aliases, append its name, description, and nicknames to the total_gene list
                total_gene.append([gene[1], gene[3], gene[2].split("|")])
            else:
                # If the gene does not have aliases, append its name, description, and None to the total_gene list
                total_gene.append([gene[1], gene[3], None])
    # Convert the total_gene list to a dictionary with gene names as keys and gene information as values
    total_gene = {total_gene[i][0]: total_gene[i][1:] for i in range(len(total_gene))}
    # Return the total_gene dictionary
    return total_gene


def create_gene_alias_table(total_gene):
    """
    Given a dictionary of genes with their aliases, creates a hash table to store the gene names and their aliases.
    """
    # Create a hash table to store the gene names and their aliases
    hash_table = {}
    for key in list(total_gene.keys()):
        hash_table[key] = key
    for key, aliases in total_gene.items():
        if aliases[1] is not None:
            for alias in aliases[1]:
                if alias not in hash_table:
                    hash_table[alias] = key
    # Return the hash table
    return hash_table


def match_gene_name(gene_name, hash_table):
    """
    Given a gene name and a hash table of genes and their aliases, checks if the gene name matches any of the aliases.
    If a match is found, returns the gene name. Otherwise, returns None.
    """
    # Check if the gene name or any of its aliases are in the hash table
    if gene_name in hash_table:
        return hash_table[gene_name]
    else:
        return None


class MarkerGene:
    # Define a class for marker genes
    all_marker_gene = {}
    # Create a dictionary to store all marker genes
    # Key: gene name, Value: gene instance
    total_gene_info = read_total_gene_in_human()
    # Read all gene information from a file
    hash_table = create_gene_alias_table(total_gene_info)

    # Create a hash table to store gene aliases

    def __init__(self, name, cell_type=None, negative=False):
        # Initialize a marker gene instance
        if cell_type is None:
            cell_type = []
        # If no cell type is provided, set it to an empty list
        self.cell_type = cell_type
        # Set the cell type of the marker gene
        self.name = name
        # Set the name of the marker gene
        self.negative = negative
        # Set whether the marker gene is negative
        MarkerGene.all_marker_gene[self.name] = self
        # Add the marker gene to the dictionary of all marker genes
        self.nicknames = MarkerGene.total_gene_info[self.name][1]
        # Get the nicknames of the marker gene
        self.description = MarkerGene.total_gene_info[self.name][0]
        # Get the description of the marker gene

    def add_cell_type(self, cell_type):
        # Add a cell type to the marker gene
        self.cell_type.append(cell_type)

    def add_nickname(self, nickname):
        # Add a nickname to the marker gene
        self.nicknames.append(nickname)
        # Add the nickname to the list of nicknames
        MarkerGene.all_marker_gene[nickname] = self


def if_marker_gene_exist(marker_gene):
    # Check if a marker gene exists
    marker_gene = str(marker_gene)
    # Convert the marker gene to a string
    search_obj = re.match("(.*)-$", marker_gene)
    # Search for a negative marker gene
    if search_obj is None:
        marker_gene = marker_gene
        negative = False
    else:
        marker_gene = search_obj.group(1)
        negative = True
    # If the marker gene is negative, set negative to True
    if marker_gene in MarkerGene.hash_table:
        return MarkerGene.hash_table[marker_gene], negative
    else:
        return None, None


def get_marker_gene(gene):
    # Get a marker gene instance
    symbol_name, negative = if_marker_gene_exist(gene)
    # Check if the marker gene exists
    if symbol_name is not None:
        if symbol_name in MarkerGene.all_marker_gene:
            if negative == MarkerGene.all_marker_gene[symbol_name].negative:
                return MarkerGene.all_marker_gene[symbol_name]
        return MarkerGene(symbol_name, negative=negative)
    else:
        return None


class Paper:
    # Define a class for papers
    paper_cell_type = {}

    # Create a dictionary to store all cell types in a paper
    # Key: paper title, Value: list of cell types

    def __init__(self, title, journal, year: int, PMID: int, method):
        # Initialize a paper instance
        self.title = title
        # Set the title of the paper
        self.journal = journal
        # Set the journal of the paper
        self.year = year
        # Set the year of the paper
        self.PMID = PMID
        # Set the PMID of the paper
        self.method = method
        # Set the method of the paper
        Paper.paper_cell_type[title] = []
        # Add the paper to the dictionary of all papers

    def toList(self):
        # Convert the paper information to a list
        return [self.title, self.journal, self.year, self.PMID, self.method]

    def find_all_relation(self):
        # Find all relations between cell types in a paper
        for cell_type in Paper.paper_cell_type[self.title]:
            cell_type.find_father_nodes()


class CellType:
    all_cell_type = {}  # A dictionary to store all instances of the CellType class

    # all_cell_type = {cell type name+paper title: cell type instance}

    def __init__(self, name, marker_genes: list, this_paper: Paper, father_node=None, source='ScSubtype'):
        # Initialize a CellType instance
        re_obj1 = re.match("(c|hM)[0-9]{2}_(.*)", name, re.IGNORECASE)  # Use regex to match the name of the cell type
        if re_obj1 is not None:
            self.name = re_obj1.group(2)  # If the name matches the regex, set the name to the second group
        else:
            self.name = name  # Otherwise, set the name to the input name
        self.source = source  # Set the source of the cell type
        self.marker_genes = []
        # Convert input marker genes to instances of marker gene class
        for i in range(len(marker_genes)):
            marker_gene = get_marker_gene(marker_genes[i])  # Get the marker gene instance
            if marker_gene is not None:
                marker_gene.add_cell_type(self)  # Add the cell type to the marker gene instance
                self.marker_genes.append(
                    marker_gene)  # Add the marker gene instance to the list of marker genes for the cell type
        self.paper = this_paper  # Set the paper of the cell type
        if this_paper is not None:
            Paper.paper_cell_type[this_paper.title].append(
                self)  # Add the cell type to the list of cell types for the paper
        self.father_node = father_node  # Set the father node of the cell type
        self.son_nodes = []  # Initialize an empty list of son nodes for the cell type
        if father_node is not None:
            father_node.son_nodes.append(self)  # Add the cell type to the list of son nodes for the father node
        if source == "ScSubtype" and this_paper is not None:
            CellType.all_cell_type[
                self.name + str(self.paper.title)] = self  # Add the cell type to the dictionary of all cell types
        else:
            CellType.all_cell_type[self.name] = self  # Add the cell type to the dictionary of all cell types
        self.nicknames = []  # Initialize an empty list of nicknames for the cell type

    def isEqual(self, cell_type_b):
        # Check if two cell types are equal
        if self.name.upper() == cell_type_b.name.upper() or self.name in cell_type_b.nicknames:
            return True
        return False

    def add_nickname(self, nickname):
        # Add a nickname to the cell type
        self.nicknames.extend(nickname)


class TypesProfileNode(CellType):
    all_profile_node = {}  # A dictionary to store all instances of the TypesProfileNode class
    # all_profile_node = {profile node name: profile node instance}
    all_own_types = []  # A list to store all cell type instances that are owned by the TypesProfileNode class

    # all_own_types = [cell type instance]

    def __init__(self, name, marker_genes: list, this_paper: Paper = None, marker_genes_weight=None, father_node=None):
        # Initialize a TypesProfileNode instance
        super().__init__(name, marker_genes, this_paper, father_node,
                         source="ProfileNode")  # Call the __init__ method of the CellType class
        if marker_genes_weight is None:
            self.marker_genes_weight = {self.marker_genes[i]: 1 / len(self.marker_genes) for i in range(
                len(self.marker_genes))}  # Set the marker gene weights for the TypesProfileNode instance
        else:
            self.marker_genes_weight = {}
            for marker_gene in marker_genes_weight:
                weight = marker_genes_weight[marker_gene]
                marker_gene = get_marker_gene(marker_gene)
                if marker_gene is not None:
                    self.marker_genes_weight[marker_gene] = weight
        if self.father_node is not None:
            self.level = self.father_node.level + 1  # Set the level of the TypesProfileNode instance
        else:
            self.level = 0
        self.son_nodes = []  # Initialize an empty list of son nodes for the TypesProfileNode instance
        self.own_types = []  # Initialize an empty list of own types for the TypesProfileNode instance
        TypesProfileNode.all_profile_node[
            self.name] = self  # Add the TypesProfileNode instance to the dictionary of all TypesProfileNode instances

    def add_marker_gene(self, marker_gene: MarkerGene, weight: float):
        # Add marker genes to the TypesProfileNode instance
        if marker_gene in self.marker_genes:
            self.marker_genes_weight[marker_gene] += weight
        else:
            self.marker_genes.append(marker_gene)
            self.marker_genes_weight[marker_gene] = weight

    # Join a cell type to the profile node.
    def join(self, cell_type: CellType):
        self.own_types.append(cell_type)
        TypesProfileNode.all_own_types.append(cell_type)


def get_source(all_profile, all_cell_type):
    for profile_node in all_profile:
        for cell_type in all_cell_type:
            if cell_type.isEqual(profile_node):
                profile_node.join(cell_type)


def get_profile_node(name, marker_genes_str: list, marker_genes_str_weight: dict, father_node=None):
    # Get the TypesProfileNode instance with the input name
    # If the TypesProfileNode instance does not exist, create a new one
    # If the TypesProfileNode instance exists, add the marker genes to the instance
    all_name = {}
    for profile_node in TypesProfileNode.all_profile_node:
        all_name[profile_node] = profile_node
        nicknames = TypesProfileNode.all_profile_node[profile_node].nicknames
        for nickname in nicknames:
            all_name[nickname] = profile_node
    if name in all_name:
        name = all_name[name]
        for marker_gene in marker_genes_str:
            weight = marker_genes_str_weight[marker_gene]
            marker_gene = get_marker_gene(marker_gene)
            if marker_gene is not None:
                TypesProfileNode.all_profile_node[name].add_marker_gene(marker_gene, weight)
        return TypesProfileNode.all_profile_node[name]
    else:
        return TypesProfileNode(name, marker_genes_str, marker_genes_weight=marker_genes_str_weight,
                                father_node=father_node)


def read_from_standard_profile():
    with open("crc_marker_profile.csv", encoding='UTF-8', newline='') as f:
        csv_reader = csv.reader(f)
        read_in = list(csv_reader)
    all_profile = []
    for i in range(0, len(read_in), 5):
        name = read_in[i][0]
        if read_in[i + 3][0] == 'None':
            father_node = None
        else:
            father_node = TypesProfileNode.all_profile_node[read_in[i + 3][0]]
        nicknames = read_in[i + 4]
        marker_genes_weight = []
        for j in range(0, len(read_in[i + 2]), 2):
            marker_genes_weight.append([read_in[i + 1][int(j / 2)], float(read_in[i + 2][j + 1])])
        # Sort marker_genes_weight by value
        marker_genes_weight = sorted(marker_genes_weight, key=lambda x: x[1], reverse=True)
        # Delete marker genes with weight less than 0.05
        marker_genes_weight = [x for x in marker_genes_weight if x[1] >= 0.05]
        marker_genes_weight = dict(marker_genes_weight)
        marker_genes = list(marker_genes_weight.keys())
        this_profile = get_profile_node(name, marker_genes, marker_genes_weight, father_node=father_node)
        if this_profile not in all_profile:
            all_profile.append(this_profile)
        this_profile.add_nickname(nicknames)
    return all_profile


def normalize_marker_gene_weight(marker_gene_weight):
    total = 0
    for key in marker_gene_weight:
        total += marker_gene_weight[key]
    for key in marker_gene_weight:
        marker_gene_weight[key] /= total
        marker_gene_weight[key] *= 100
    return marker_gene_weight


def rank_marker_gene(marker_gene_weight, coefficient):
    # Sort marker_genes_weight by value
    marker_gene_weight = sorted(marker_gene_weight.items(), key=lambda x: x[1], reverse=True)
    marker_gene_weight = dict(marker_gene_weight)
    for i in range(len(marker_gene_weight)):
        marker_gene_weight[list(marker_gene_weight.keys())[i]] = coefficient ** i
    return marker_gene_weight


def normalize_self_marker_gene(marker_genes_weight):
    scaled_marker_gene_weight = {}
    # Sort marker_genes_weight by value
    marker_genes_weight = sorted(marker_genes_weight.items(), key=lambda x: x[1], reverse=True)
    avg = sum(list(map(lambda x: x[1], marker_genes_weight))) / len(marker_genes_weight)
    std = np.std(list(map(lambda x: x[1], marker_genes_weight)))
    # Scale marker_genes_weight
    for i in range(len(marker_genes_weight)):
        scaled_marker_gene_weight[marker_genes_weight[i][0]] = (marker_genes_weight[i][1] - avg) / std
    return scaled_marker_gene_weight


def scale_among_same_marker_gene(marker_gene_list):
    for i in range(len(marker_gene_list)):
        total = 0
        for j in range(len(marker_gene_list[i])):
            if marker_gene_list[i][j] > 0:
                total += 1
        if total > 1:
            for j in range(len(marker_gene_list[i])):
                marker_gene_list[i][j] /= total
    return marker_gene_list


def enrich(this_cell_type: TypesProfileNode, coefficient=0.1):
    next_cell_type = this_cell_type.son_nodes
    if len(next_cell_type) != 0:
        for child_cell_type in next_cell_type:
            enrich(child_cell_type, coefficient)
            for marker_gene, weight in child_cell_type.marker_genes_weight.items():
                this_cell_type.add_marker_gene(marker_gene, weight * coefficient)


def marker_gene_matrix(sub_cell_type, coefficient):
    rowName = []
    # row name is marker gene name
    colName = []
    # column name is cell type
    for cell_type in sub_cell_type:
        colName.append(cell_type)
        for marker_gene in list(map(lambda x: x.name, cell_type.marker_genes)):
            if marker_gene not in rowName:
                rowName.append(marker_gene)
    marker_gene_list = [[0 for i in range(len(colName))] for j in range(len(rowName))]
    for cell_type in sub_cell_type:
        col = colName.index(cell_type)
        for i in range(len(cell_type.marker_genes)):
            marker_gene = cell_type.marker_genes[i]
            row = rowName.index(marker_gene.name)
            marker_genes_weight = rank_marker_gene(cell_type.marker_genes_weight, coefficient)
            if marker_gene.negative:
                marker_gene_list[row][col] = -(marker_genes_weight[marker_gene])
            else:
                marker_gene_list[row][col] = marker_genes_weight[marker_gene]
    marker_gene_list = scale_among_same_marker_gene(marker_gene_list)
    return marker_gene_list, rowName, colName


def most_similar(sub_cell_type, query_marker_gene_list, score_cutoff, rank_coefficient):
    marker_gene_list, rowName, colName = marker_gene_matrix(sub_cell_type, rank_coefficient)
    query_list = [1 for i in range(len(query_marker_gene_list))]
    reference_list = []
    for i in range(len(query_marker_gene_list)):
        if query_marker_gene_list[i][0].name in rowName:
            reference_list.append(marker_gene_list[rowName.index(query_marker_gene_list[i][0].name)])
        else:
            reference_list.append([0 for i in range(len(colName))])
    query_matrix = np.array([query_list])
    result_matrix = np.dot(query_matrix, np.array(reference_list))
    max_match = result_matrix.argmax()
    print(list(map(lambda x: x.name, colName)))
    print(result_matrix)
    if result_matrix.max() <= score_cutoff:
        return "Unaligned cell"
    else:
        return colName[max_match]


def search_through_level(query_marker_gene_list, score_cut_off, depth_cut_off, rank_coefficient, start_node):
    profile_in_this_level = [TypesProfileNode.all_profile_node[start_node]]
    depth = 0
    while len(profile_in_this_level) != 0:
        profile_in_next_level = []
        for profile in profile_in_this_level:
            if not profile.son_nodes:
                if len(profile_in_this_level) == 1:
                    return profile.name
            else:
                depth += 1
                for child in profile.son_nodes:
                    if child not in TypesProfileNode.all_own_types:
                        profile_in_next_level.append(child)
        if len(profile_in_next_level) <= 0:
            return "Unaligned cell"
        else:
            most_similar_cell_type = most_similar(profile_in_next_level, query_marker_gene_list, score_cut_off,
                                                  rank_coefficient)
            if depth >= depth_cut_off:
                return most_similar_cell_type.name
            if most_similar_cell_type == "Unaligned cell":
                profile_in_this_level = profile_in_next_level
            else:
                profile_in_this_level = [most_similar_cell_type]


def read_query(file, all_cell_type, rank_coefficient=0.9, score_cut_off=0, depth_cut_off=2,
               enrichment_coefficient=0.5, start_node='Unaligned cell'):
    enrich(all_cell_type[0], enrichment_coefficient)
    # Read the query file and create a list of gene of cell clusters
    with open(file, newline='') as f:
        reader = csv.reader(f)
        query = list(reader)
    query_cluster = []
    for i in range(1, len(query)):
        cluster = int(query[i][5])
        lfc = float(query[i][1])
        gene = query[i][6]
        marker_gene = get_marker_gene(gene)
        if marker_gene is not None:
            marker_gene = [marker_gene, lfc]
            if cluster != len(query_cluster) - 1:
                query_cluster.append([marker_gene])
            else:
                query_cluster[-1].append(marker_gene)
    # Call the `most_similar` function for each gene cluster and append the result to a list
    result = []
    for i in range(len(query_cluster)):
        result.append(search_through_level(query_cluster[i], score_cut_off, depth_cut_off, rank_coefficient, start_node))
    # Return the list of cell types that match the input gene clusters
    return result


def output(result, file):
    with open(file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(result)


def read_from_standard():
    with open("crc_marker.csv", encoding='UTF-8', newline='') as f:
        csv_reader = csv.reader(f)
        input = list(csv_reader)
    this_line = 0
    total_cell_types = []
    while this_line < len(input):
        if str(input[this_line][0] + input[this_line + 2][0]) in CellType.all_cell_type:
            this_line += 5
            continue
        title = input[this_line + 2][0]
        journal = input[this_line + 2][1]
        year = int(input[this_line + 2][2])
        PMID = int(input[this_line + 2][3])
        method = input[this_line + 2][4]
        paper = Paper(title, journal, year, PMID, method)
        father_node = input[this_line + 3]
        if father_node[0] != 'None':
            father_node = CellType.all_cell_type[str(father_node[0]) + str(father_node[1])]
            total_cell_types.append(CellType(input[this_line][0], input[this_line + 1], paper, father_node))
        else:
            total_cell_types.append(CellType(input[this_line][0], input[this_line + 1], paper))
        this_line += 5
    return total_cell_types


def main(avg):
    profile = read_from_standard_profile()
    read = read_from_standard()
    get_source(profile, read)
    result = read_query(avg[0], profile, avg[2], avg[3], avg[4], avg[5], avg[6])
    print(result)
    output(result, avg[1])


if __name__ == "__main__":
    avg = (r"\\10.20.5.147\data417\dingcy\rstudio\2020 Lineage-dependent gene expression programs influence the "
           r"immune landscape of colorectal cancer\crc_CD4_cell_markers_resolution=2.csv",
           r"\\10.20.5.147\data417\dingcy\rstudio\2020 Lineage-dependent gene expression programs influence the "
           r"immune landscape of colorectal cancer\crc_CD4_cell_markers_result.csv", 0.9, 0, 10, 0.5, 'CD4 T cell'
           )
    main(avg)
