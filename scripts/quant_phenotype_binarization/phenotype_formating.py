import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# CONSTANT VARIABLES; adjust for your specific phenotype / genotype files characteristics
phenotype_file = 'aggregated_mean.csv'
phenotype_file_header_present = 'Y'
phenotype_file_sampleID_column_index = 0 # index of accesion numbers column in phenotype file (indexing starts from 0)
phenotype_file_values_column_index = 1 # index of phenotype values column in phenotype file (indexing starts from 0)

genotype_file = 'filtered_data.csv'
genotype_file_header_present = 'Y'
genotype_file_sampleID__column_index = 1 # index of accesion numbers column in genotype file (indexing starts from 0)
genotype_file_values_column_index = 2 # index of binarised genotype values column in genotype file (indexing starts from 0)
genotype_file_type_column_index = 3 # index of type column in genotype file (landrace, cultivar ...)

separator = ',' # both genotype and phenotype file should use same separator

function_list = []
info_strings = []
def register_with_info(info):
    ''''Decorator factory that creates a decorator which registers a function and an info string. '''
    def decorator(function):
        function_list.append(function)
        info_strings.append(info)
        return function
    return decorator

# FUNCTIONS

def read_data(file, indices, column_names, separator, header):
    """
    Read specific columns from a CSV file and assign new column names.

    Parameters:
    - file (str): genotype_file / Phenotype_file
    - indices (list): List of required columns indices to extract.
    - column_names (list): List of new names to assign to the columns.
    - separator (str): Delimiter used in file.
    - header (str): 'Y' if the file includes a header row, otherwise 'N'.

    Returns:
    - A pandas DataFrame with selected columns and specified column names.
    """
    header = 0 if header == 'Y' else None
    df = pd.read_csv(file, sep=separator, header=header, usecols=indices)
    df.columns = column_names
    return df

def read_phenotype_data(file, sample_index, phenotype_index, separator, header):
    '''Read and process the phenotype file data, extracting the required columns into a pandas DataFrame.'''
    df = read_data(file, [sample_index, phenotype_index], ['Sample', 'Phenotype'], separator, header)
    df['Phenotype'] = pd.to_numeric(df['Phenotype'], errors='coerce')  # Convert to numeric, coercing errors to NaN
    return df

def read_genotype_data(file, sample_index, genotype_index, type_index, separator, header):
    '''Read and process the genotype file data, extracting the required columns into a pandas DataFrame.'''
    return read_data(file, [sample_index, genotype_index, type_index], ['Sample', 'Genotype', 'type'], separator, header)

# Data statistics
def compute_statistics(df):
    '''Compute statistical parameters for phenotype data and store them into dictionary'''
    min_value = df.Phenotype.min()
    max_value = df.Phenotype.max()
    range_value = max_value - min_value

    return {
        'min': min_value,
        'max': max_value,
        'mean': df.Phenotype.mean(),
        'median': df.Phenotype.median(),
        'first_quartile': df.Phenotype.quantile(0.25),
        'third_quartile': df.Phenotype.quantile(0.75),
        'range': range_value
    }

# quantitative phenotype binarisation:

# 1st approach to binarisation
@register_with_info("Average accuracy using median value for binarisation: ")
def binarisation_by_median(df, statistics):
    '''Binarize the phenotype data based on the median.'''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['median']).astype(int) # overwrite Phenotype values 
    return df_copy

# 2nd approach to binarisation
@register_with_info("Average accuracy using binarisation by extremes defined by 1st quartile and 3rd quartile: ")
def binarisation_by_extremes(df, statistics):
    '''
    Binarize the phenotype data based on defined extreme intervals:
    - First extreme interval: values from the minimum to the first quartile are classified as 0.
    - Second extreme interval: values from the third quartile to the maximum are classified as 1.
    Values outside these intervals are set to -2.
    '''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original

    conditions = [
    df_copy['Phenotype'] < statistics['first_quartile'], # zeros
    df_copy['Phenotype'] > statistics['third_quartile'] # ones
    ]
    choices = [0, 1] # Values to assign based on the conditions

    # Apply conditions to the data, setting out-of-range values to -2
    df_copy['Phenotype'] = np.select(conditions, choices, default = -2)
    return df_copy

# 3rd approach to binarisation
@register_with_info(f"Average accuracy using 1st quartile value for binarisation: ")
def binarisation_by_first_quartile(df, statistics):
    '''Binarize the phenotype data based on the first quartile.'''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['first_quartile']).astype(int)
    return df_copy

# 4th approach to binarisation
@register_with_info("Average accuracy using 3rd quartile value for binarisation: ")
def binarisation_by_third_quartile(df, statistics):
    '''Binarize the phenotype data based on the third quartile.'''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['third_quartile']).astype(int)
    return df_copy

# 5th approach
def binarisation_combinations (genotype_df, phenotype_df):
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
    for segment_zeros in range(0, 100, interval):
        for segment_ones in range(0, 100, interval):
            if segment_zeros == segment_ones:
                continue
            else:
                zeros_quantil_start = segment_zeros / 100
                zeros_quantil_end = (segment_zeros + interval) / 100
                ones_quantil_start = segment_ones / 100
                ones_quantil_end = (segment_ones + interval) / 100

                conditions = [
                    (df_copy['Phenotype'] > df_copy.Phenotype.quantile(zeros_quantil_start)) & (df_copy['Phenotype'] < df_copy.Phenotype.quantile(zeros_quantil_end)),
                    (df_copy['Phenotype'] >  df_copy.Phenotype.quantile(ones_quantil_start)) & (df_copy['Phenotype'] < df_copy.Phenotype.quantile(ones_quantil_end))
                    ]
                choices = [0, 1]  # Corresponding values for the conditions above
                
                df_copy2 = df_copy.copy() # Create second copy to store binarised phenotype. 1st copy is used for other combinations
                df_copy2['Phenotype'] = np.select(conditions, choices, default = -2)

                #df_copy['Binarized Phenotype'] = np.select(conditions, choices, default=-2)
                # Forming the string
                info_string = f"Borders of 0s: {zeros_quantil_start}, {zeros_quantil_end} \nborders of 1s: {ones_quantil_start}, {ones_quantil_end}"

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

def average_accuracy_from_dataframes(genotype_dataframe, phenotype_dataframe):
    '''
    Compute the average accuracy of genotype-phenotype associations for mutant and wild type classes.
    accuracy formula: average accuracy = (MUT_acc + WT_acc) / 2
    MUT_acc - accuracy of correct mutant genotye - mutant fenotype association (1/1)
    WT_acc accuracy of correct wild type genotype - wild type phenotype association (0/0)
    '''
    # Merge the dataframes on the 'GRIN_Accession' column
    merged_df = pd.merge(genotype_dataframe, phenotype_dataframe, on='Sample')

    # formulas needed for average accuracy
    N_WTcorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 0)])
    N_WTincorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 1)])
    acc_WT = (N_WTcorrect / (N_WTcorrect + N_WTincorrect)) * 100

    N_MUTcorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 1)])
    N_MUTincorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 0)])
    acc_MUT = (N_MUTcorrect / (N_MUTcorrect + N_MUTincorrect)) * 100

    # average accuracy formula
    return (acc_MUT + acc_WT) / 2

def histogram(df):
    ''' Create and plot the histogram for the phenotype data '''
    bin_edges = np.histogram_bin_edges(phenotype_df['Phenotype'].dropna(), bins='auto')  # calculates the optimal number of bins based on the data

    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(phenotype_df['Phenotype'], bins=bin_edges, color='#93BFCF', edgecolor='#6096B4', histtype='stepfilled', linewidth=1.5, alpha=0.8)

    # Add grid, labels, and title
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Phenotype Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of trait values')

    # Save and show the plot
    plt.savefig("histogram.png", dpi=300, bbox_inches='tight')
    plt.show()


def boxplot(genotype_df, phenotype_df):
    '''
    Visualize the distribution of phenotype values for accessions with 
    reference and alternative genotypes.
    '''
    genotype_df.to_csv('genotype.csv',  index=False)
    phenotype_df.to_csv('phenotype.csv',  index=False)

    merged_df = pd.merge(genotype_df, phenotype_df, on='Sample')
    merged_df.to_csv('boxplot_data.csv',  index=False)


    # Create a figure and a set of subplots
    fig, ax = plt.subplots()

    # Predefining required properties of boxplot (the box, median line, whiskers and caps appearance)
    box_properties = dict(linestyle='-', linewidth=3, color='#6096B4', facecolor='#93BFCF', alpha=0.4)
    median_properties = dict(linestyle='-', linewidth=2.5, color='firebrick')
    whisker_properties = dict(linestyle='-', linewidth=3, color='#6097B466')
    caps_properties = dict(linestyle='-', linewidth=3, color='#6097B466')

    # This creates the boxplot with the properties for box and median
    # 'patch_artist' must be set on True for properties like facecolor and alpha to work
    bp = merged_df.boxplot(column='Phenotype', by='Genotype', ax=ax, boxprops=box_properties, medianprops=median_properties, whiskerprops = whisker_properties, capprops = caps_properties, showfliers=False, patch_artist=True)

    conditions = [
                (merged_df['type'] == 'Elite'),
                (merged_df['type'] == 'G. soja'),
                (merged_df['type'] == 'Genetic'),
                (merged_df['type'] == 'Landrace')
                ]
    choices = ['blue', 'green', 'red', 'purple']  # Corresponding values for the conditions above
    
    merged_df['color'] = np.select(conditions, choices, default = 'black')

    # Adding jitter and  generating correct boxplots columns labels
    # create 2 groups based on genotype 
    # - Group 1: reference genotype rows
    # - Group 2: alternative genotype rows
    # enumerate encapsulates all objects (DataFrameGroupBy) merged_df.groupby returns (In this case 2 objects: Both objects contain two values (key, group)
    name_setter = [] # list for labels describing a boxplot ('reference', 'alternative')
    name_setter_index = [] # correct indexes for titles, boxplots are indexed from 1, because index 0 belongs to y ax
    for i, (key, group) in enumerate(merged_df.groupby('Genotype')):
         # Append labels and indices
        if key == 0:
            name_setter.append(f'Reference Genotype ({key})')
            name_setter_index.append(i + 1)
        elif key == 1:
            name_setter.append(f'Alternative Genotype ({key})')
            name_setter_index.append(i + 1)

        # jitered x values:
        x = np.random.normal(i + 1, 0.017, size=len(group))
        ax.scatter(x, group['Phenotype'], alpha=0.5, color=group['color'], edgecolor='none')

    # Additional plot formatting
    ax.set_title('Phenotype Distribution for accesions with reference and alternative genotype')
    plt.suptitle('')  # Removes the 'Group By' title
    plt.xlabel('')
    plt.ylabel('Phenotype Value')
    plt.xticks(name_setter_index, name_setter)  # Sets custom x-axis labels

    # Show the plot
    plt.show()

# _____________________________________________ MAIN ___________________________________________

# load the data from phenotype and genotype file into dataframe

phenotype_df = read_phenotype_data(phenotype_file, phenotype_file_sampleID_column_index, phenotype_file_values_column_index, separator, phenotype_file_header_present)
binarised_genotype_df = read_genotype_data(genotype_file, genotype_file_sampleID__column_index, genotype_file_values_column_index, genotype_file_type_column_index, separator, genotype_file_header_present)


stats = compute_statistics(phenotype_df)

for function, info in zip(function_list, info_strings):
    binarised_phenotype_df = function(phenotype_df, stats)
    avg = average_accuracy_from_dataframes(binarised_genotype_df, binarised_phenotype_df)
    print(f'{info}{avg:.2f}%')

binarisation_combinations(binarised_genotype_df, phenotype_df)
# OUTPUT
histogram(phenotype_df)
boxplot(binarised_genotype_df, phenotype_df)
print(f"The phenotype file contains {len(phenotype_df.Phenotype)} samples.")
print(f"The MINIMAL VALUE of the trait is:   {stats['min']:.2f}.")
print(f"The MAXIMAL VALUE of the trait is:   {stats['max']:.2f}.")
print(f"The AVERAGE VALUE of the trait is:   {stats['mean']:.2f}.")
print(f"The MEDIAN VALUE of the trait is:    {stats['median']:.2f}.")
print(f"The 1ST QUARTILE of the trait is:    {stats['first_quartile']:.2f}.")
print(f"The 3RD QUARTILE of the trait is:    {stats['third_quartile']:.2f}.")