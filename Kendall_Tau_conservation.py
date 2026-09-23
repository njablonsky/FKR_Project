import pandas as pd
from scipy.stats import kendalltau
from itertools import combinations

# 1. Load the data
df = pd.read_csv("marker_location.txt", sep="\t")

# 2. Extract ordered subfamily sequences for each genome
genome_sequences = {}
for genome, group in df.groupby('Genome'):
    # Sort by the physical Location to get the true sequence order
    ordered_genes = group.sort_values('Location')['Subfam'].tolist()
    genome_sequences[genome] = ordered_genes

genomes = sorted(list(genome_sequences.keys()))
results = []

# 3. Iterate through every pair of genomes
for g1, g2 in combinations(genomes, 2):
    seq1_full = genome_sequences[g1]
    seq2_full = genome_sequences[g2]
    
    # Get the shared markers (order from seq1)
    set2 = set(seq2_full)
    common_genes = [gene for gene in seq1_full if gene in set2]
    
    n = len(common_genes)
    if n < 3:
        # P-value calculation is only robust with 3 or more points
        continue
        
    # Get order of shared markers in seq2
    set1 = set(seq1_full)
    common_in_seq2 = [gene for gene in seq2_full if gene in set1]
    
    # Map common genes to their ranks in the query sequence
    rank_map = {gene: i for i, gene in enumerate(common_genes)}
    ranks_in_seq2 = [rank_map[gene] for gene in common_in_seq2]
    
    # 4. Calculate Kendall Tau and p-value
    tau, p_val = kendalltau(range(n), ranks_in_seq2)
    
    results.append({
        'Genome_A': g1,
        'Genome_B': g2,
        'Shared_Markers': n,
        'Kendall_Tau': round(tau, 4),
        'P_Value': p_val
    })

# 5. Save results to a CSV file
output_df = pd.DataFrame(results)
output_df.to_csv("pairwise_conservation_stats.csv", index=False)
print("Analysis complete. Results saved to pairwise_conservation_stats.csv")
