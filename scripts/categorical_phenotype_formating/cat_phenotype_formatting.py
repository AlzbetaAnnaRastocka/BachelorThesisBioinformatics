import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_data(file, indices, column_names, separator, header = -1):
    """Read specific columns from a CSV file and assign new column names."""
    header = header if header >= 0 else None
    df = pd.read_csv(file, sep=separator, header=header, usecols=indices)
    df.columns = column_names
    return df

def read_phenotype_data(file, sample_index, phenotype_index, separator, header, phenotype_categories):
    """Read and process the phenotype file data, extracting the required columns into a pandas DataFrame."""

    df = read_data(file, sample_index + phenotype_index, ['accession_prefix', 'accession_number', 'accession_suffix','Phenotype'], separator, header)

    df['accession_prefix'] = df['accession_prefix'].fillna('')
    df['accession_number'] = df['accession_number'].fillna('')
    df['accession_suffix'] = df['accession_suffix'].fillna('')

    df['Sample'] = df['accession_prefix'].astype(str) + df['accession_number'].astype(str) + df['accession_suffix'].astype(str)

    new_df = df[['Sample', 'Phenotype']]
    new_df = new_df.drop_duplicates()
    new_df = new_df.dropna(subset=['Phenotype'])

    # Define the categorical type with the specific order
    phenotype_cats = pd.CategoricalDtype(phenotype_categories, ordered=True)
    new_df['Phenotype'] = new_df['Phenotype'].astype(phenotype_cats)  # Convert column to categorical type

    return new_df

def read_genotype_data(file, sample_index_GRIN, binarised_genotype_index, separator, header):
    """Read and process the phenotype file data, extracting the required columns into a pandas DataFrame."""

    df = read_data(file, sample_index_GRIN + binarised_genotype_index, ['Sample', 'Genotype'], separator, header)
    return df

def histogram(df):
    """Create and plot the histogram for the phenotype data."""
    plt.figure(figsize=(10, 6))
    # Plot directly from the categorical data, pandas will handle the categories correctly
    ax = df['Phenotype'].value_counts().sort_index().plot(kind='bar', color='#93BFCF', edgecolor='#6096B4', linewidth=1.5, alpha=0.5)
    plt.xlabel('Phenotype Groups')
    plt.ylabel('Frequency')
    plt.title('Distribution of Phenotype Groups')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xticks(rotation=0)  # Change rotation parameter to 0 for horizontal labels

    for p in ax.patches:  # ax.patches is a list of rectangles, with each rectangle representing a bar
        ax.annotate(f'{int(p.get_height())}',  # The text to display
                     (p.get_x() + p.get_width() / 2., p.get_height()),  # Position (x, y)
                     ha='center',  # Horizontal alignment (centering the text)
                     va='center',  # Vertical alignment
                     xytext=(0, 9),  # Distance from the top of the bar
                     textcoords='offset points')  # How the text is positioned
    plt.show()

# Phenotype binarisation
def phenotype_binarisation(df):
    df_copy = df.copy()
    conditions = [
        df_copy['Phenotype'].isin(['IV']),
        df_copy['Phenotype'].isin(['000', '00', '0']),
    ]
    choices = [0, 1]
    df_copy['Phenotype'] = np.select(conditions, choices, default = -2)

    return df_copy
    

def average_accuracy_from_dataframes(genotype_dataframe, phenotype_dataframe):
    '''
    Compute the average accuracy of genotype-phenotype associations for mutant and wild type classes.
    accuracy formula: average accuracy = (MUT_acc + WT_acc) / 2
    MUT_acc - accuracy of correct mutant genotye - mutant fenotype association (1/1)
    WT_acc accuracy of correct wild type genotype - wild type phenotype association (0/0)
    '''
    # Merge the dataframes on the 'GRIN_Accession' column
    merged_df = pd.merge(genotype_dataframe, phenotype_dataframe, on='Sample')
    merged_df.to_csv('merged.csv', index=False)

    # formulas needed for average accuracy
    N_WTcorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 0)])
    N_WTincorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 1)])
    acc_WT = 0 if (N_WTcorrect + N_WTincorrect) == 0 else (N_WTcorrect / (N_WTcorrect + N_WTincorrect)) * 100

    N_MUTcorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 1)])
    N_MUTincorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 0)])
    acc_MUT = 0 if (N_MUTcorrect + N_MUTincorrect) == 0 else (N_MUTcorrect / (N_MUTcorrect + N_MUTincorrect)) * 100

    # average accuracy formula
    return (acc_MUT + acc_WT) / 2



def binarisation_combinations (genotype_df, phenotype_df, categories):
    '''
    Evaluates and identifies the best binarization method for phenotype data based on quantile intervals to maximize the correspondence with binarized genotype data.
    The function explores all possible pairs of quantile intervals (excluding identical pairs) to determine which combination yields the highest accuracy. 
       The size of the interval for quantile division is determined by the 'interval' variable, which accepts values that evenly divide 100 (e.g., 5, 10, 20, 25, 50).
       We store the average accuracy values of all the combinations in a dictionary and then pick the best value.   
    '''
    df_copy = phenotype_df.copy() # df_copy is a copy of the phenotype DataFrame to avoid modifying the original
    combinations_for_binarisation_dict = dict()
    interval = 10
    dictionary_index = 0
    for segment_zeros in categories:
        for segment_ones in categories:
            if segment_zeros == segment_ones:
                continue
            else:

                conditions = [
                    (df_copy['Phenotype'] == segment_zeros),
                    (df_copy['Phenotype'] == segment_ones)
                    ]
                choices = [0, 1]  # Corresponding values for the conditions above
                
                df_copy2 = df_copy.copy() # Create second copy to store binarised phenotype. 1st copy is used for other combinations
                df_copy2['Phenotype'] = np.select(conditions, choices, default = -2)

                #df_copy['Binarized Phenotype'] = np.select(conditions, choices, default=-2)
                # Forming the string
                info_string = f"Category for 0s: {segment_zeros}, \nCategory for 1s: {segment_ones}"

                # Updating the dictionary with a tuple containing the dataframe slice and the string
                combinations_for_binarisation_dict.update({ dictionary_index: (df_copy2, info_string)})
                dictionary_index += 1
                
    accuracies = dict()
    for key, value in combinations_for_binarisation_dict.items():
        # value is a tuple where value[0] is the dataframe and value[1] is the string
        accuracy = average_accuracy_from_dataframes(genotype_df, value[0])
        accuracies.update({key: accuracy})

    key_with_highest_value = max(accuracies, key=accuracies.get)
    value_best = combinations_for_binarisation_dict[key_with_highest_value]
    print(f'Key with the highest value: {key_with_highest_value}, Value: {accuracies[key_with_highest_value]}, {value_best[1]}')




# ______________________________ MAIN ______________________________s

# Example usage
phenotype_file = 'grin.csv'
phenotype_file_header = 2
phenotype_file_sample_index = [1,2,3]  # index of accession numbers column
phenotype_file_phenotype_index = [4]  # index of phenotype values column
separator = ','
phenotype_categories = ['000', '00', '0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']
number_of_categories = len(phenotype_categories)

binarised_genotype_file = "GRIN_genotype.csv"
genotype_file_header = 0

# Read the data
phenotype_df = read_phenotype_data(phenotype_file, phenotype_file_sample_index, phenotype_file_phenotype_index, separator, phenotype_file_header, phenotype_categories)
binarised_phenotype_df = phenotype_binarisation(phenotype_df)

binarised_genotype_df = read_genotype_data(binarised_genotype_file, [0], [1], ',', genotype_file_header)


print(average_accuracy_from_dataframes(binarised_genotype_df, binarised_phenotype_df))

binarisation_combinations(binarised_genotype_df, phenotype_df, phenotype_categories)

# Plot the histogram
histogram(phenotype_df)
