import os

def read_file(directory): # Opens all files and returns the gene bases upstream of start codon 
    upstream = {}
    for file_name in os.listdir(directory):
        if file_name.endswith(".txt"):  
            gene_name = file_name.split('.')[0]  # Uses file name as gene name
            with open(os.path.join(directory, file_name), 'r') as file:
                gene = file.read().replace(" ", "")  # Removes spaces
                upstream[gene_name] = gene[:17]  # Stores the first 17 bases upstream of the start codon 
    return upstream

def count_mis(shine):
    seq = "AGGAGG" # The sequence we're looking for 
    lowest_i = []  # The index at which the lowest mismatch sequence starts 
    lowest = float('inf')  # Start with an infinitely high mismatch count as default
    best_shine = None #Start with None as default 

    for i in range(len(shine) - len(seq) + 1):  # Iterates through each base 
        mis_count = 0
        for j in range(len(seq)):  # Iterates through 6 bases to see number of mismatches 
            if shine[i + j] != seq[j]:
                mis_count += 1
        
        if mis_count < lowest: 
            lowest = mis_count
            lowest_i = [i]  # Resets list with the new best index
            best_shine = shine[i:i + len(seq)]  # Stores best matching sequence
        elif mis_count == lowest:
            lowest_i.append(i)  #If the nb of mismatches is equal to lowest, add the result to the list 

    # Finds the separation distance 
    if len(lowest_i) == 1:
        separation = len(shine) - lowest_i[0] - len(seq)  # Separation to start codon
    else:
        best_shine = None  # If there are multiple matches, set shine sequence to None
        separation = None  # If there are multiple matches, set separation sequence to None

    return best_shine, lowest, separation

def analyze_genes(directory):
    
    gene_data = read_file(directory) # Open files and return 17 upstream bases of start codon 
    
    # The header for the "table"
    print(f"{'Gene':<10} {'Shine':<10} {'Mismatches':<12} {'Separation':<10} {'17 upstream bases':<20}")
    
    #Iterate through the gene data retrieved 
    for gene_name, upstream_seq in gene_data.items():
        
        #Store the results of count_mis()
        best_shine, mismatches, separation = count_mis(upstream_seq)
        
        # Turn best_shine and separation into readable strings 
        if best_shine:
            shine = best_shine
        else:
            shine = "None"

        if separation != None:
            sep = separation
        else:
            sep = "None"
        
        # Add rows to the table with the results 
        print(f"{gene_name:<10} {shine:<10} {mismatches:<12} {sep:<10} {upstream_seq:<20}")


# Call the function to analyze and print results, give the directory of the gene files 
analyze_genes("C:/Users/Cindy/OneDrive/Desktop/CPS 501 (BIO)/a3")
