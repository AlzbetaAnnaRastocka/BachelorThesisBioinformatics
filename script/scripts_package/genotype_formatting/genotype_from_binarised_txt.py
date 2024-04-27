import pandas as pd

def main():
    genotype = 'genotype.txt'
    genotype_with_GRIN = 'Soybean_accession_data.csv'

    column_names = ['Accession', 'Genotype']
    genotype_df = pd.read_csv(genotype, delimiter='\t', names = column_names)  # Adjust this based on the file delimiter
    genotype_df = genotype_df.dropna(subset=['Genotype'])
    genotype_df['Genotype'] = genotype_df['Genotype'].astype(int)

    columns_to_load = ['Accession', 'GRIN_Accession', 'Classification']
    genotype_GRIN_accessions_df = pd.read_csv(genotype_with_GRIN, usecols=columns_to_load)
    genotype_GRIN_accessions_df.rename(columns={'GRIN_Accession': 'Sample', 'Classification' : 'type'}, inplace=True)

    merged_df = pd.merge(genotype_GRIN_accessions_df, genotype_df, on='Accession')
    new_df = merged_df[['Sample', 'Genotype']]
    new_df.to_csv("genotype_binarised.csv", index=False)