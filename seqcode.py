import os

#Opens the gene file and returns it as a string 
def open_gene(filepath):
    sequence = ""
    with open(filepath, 'r') as f:
        #Seperates the gene 
        for line in f:
            sequence += line.strip()  
    return sequence

#Opens the variant file and opens it as a string 
def open_fasta(filepath):
    sequence = ""
    with open(filepath, 'r') as f:
        for line in f:
            #Seperates the gene 
            if not line.startswith('>'):  
                sequence += line.strip()  
    return sequence

#Calculates the allowed number of mismatches 
def allowed_mismatches(gene_length):
    allowed = int(0.01 * gene_length)
    return max(5, min(allowed, 20))

#Coompares one gene and variant sequence
def compare_sequences(reference_gene, variant_seq, gene_name):
    gene_len = len(reference_gene)
    allowed = allowed_mismatches(gene_len)
    mutation_list = []

    found_match = False
    for i in range(len(variant_seq) - gene_len + 1):
        mismatches = 0
        mutations = []

        for j in range(gene_len):
            ref_nuc = reference_gene[j]
            var_nuc = variant_seq[i + j]
            
            if ref_nuc != var_nuc:
                mismatches += 1
                mutations.append(f"{j+1}{ref_nuc}/{var_nuc}")
            
            # If mismatches are greater than allowed, stop comparing 
            if mismatches > allowed:
                break
        
        # If within allowed mismatches, record the mutation information
        if mismatches <= allowed:
            found_match = True
            mutation_list.append({
                'start_position': i,
                'mismatch_count': mismatches,
                'mutations': mutations
            })
    
    if mutation_list:
        for mutation in mutation_list:
            print(f"\nGene: {gene_name}")
            print(f"Start position in variant: {mutation['start_position']}")
            print(f"Mismatch count: {mutation['mismatch_count']}")
            if mutation['mutations']:
                print(f"Mutations: {', '.join(mutation['mutations'])}")
            else:
                print(f"Mutations: None")
    elif not found_match:
        print(f"No matching sequence found for gene {gene_name} in the variant.")
    
    return mutation_list

#Compare all genes to a variant
def compare_genes_to_variant(ref_dir, variant_file):
    variant_seq = open_fasta(variant_file)  # Open the full variant 
    mutation_report = {}

    # Iterate through all gene files 
    for ref_gene_file in os.listdir(ref_dir):
        if ref_gene_file.endswith(".txt"):
            gene_name = ref_gene_file.split(".")[0]  # Extract gene name
            ref_gene_seq = open_gene(os.path.join(ref_dir, ref_gene_file))  # Open reference gene sequence
            mutations = compare_sequences(ref_gene_seq, variant_seq, gene_name)
            mutation_report[gene_name] = mutations
    
    return mutation_report

# Define file paths 
reference_genes_directory = "reference_genes_directory/"  # Where the genes are stored 
variant_files = [
    "PQ471628.1.fasta", "OZ195338.1.fasta", "PQ481355.1.fasta", "PQ479843.1.fasta", "OZ196968.1.fasta"
]  

# Compare each variant sequence with all reference genes
for variant_file in variant_files:
    print(f"\nComparing variant: {variant_file}")
    mutation_report = compare_genes_to_variant(reference_genes_directory, variant_file)

