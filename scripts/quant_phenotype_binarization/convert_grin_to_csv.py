import pandas as pd

path_to_file = 'grin.csv'

# formatting GRIN real phenotype data
df = pd.read_csv(path_to_file, sep = ',', skiprows=2)

df['accession_prefix'] = df['accession_prefix'].fillna('')
df['accession_number'] = df['accession_number'].fillna('')
df['accession_suffix'] = df['accession_suffix'].fillna('')

# unite 'accession_prefix' and 'accesion_number' into one 'GRIN_Accession' column to obtain single value without separator example: PI587848
df['GRIN_Accession'] = df['accession_prefix'].astype(str) + df['accession_number'].astype(str) + df['accession_suffix'].astype(str)

# grouping lines by the GRIN_Accessions and collect observed into lists for each GRIN_accession, make mean value from lists
df['mean_observation_value'] = df.groupby('GRIN_Accession')['observation_value'].transform(lambda x: pd.Series(x).mean())

df = df.drop_duplicates(subset=['GRIN_Accession', 'mean_observation_value'])
df_aggregated_mean = df[['GRIN_Accession', 'mean_observation_value']]
df_aggregated_mean.to_csv('aggregated_mean.csv', index=False)
