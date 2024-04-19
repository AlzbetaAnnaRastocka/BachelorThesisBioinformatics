import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
# _____________________________________________ read data  _____________________________________________

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

# _____________________________________________ read phenotype  _____________________________________________

def read_phenotype_data(file, sample_index, phenotype_index, separator, header):
    '''Read and process the phenotype file data, extracting the required columns into a pandas DataFrame.'''
    df = read_data(file, [sample_index, phenotype_index], ['Sample', 'Phenotype'], separator, header)
    df['Phenotype'] = pd.to_numeric(df['Phenotype'], errors='coerce')  # Convert to numeric, coercing errors to NaN
    return df

# _____________________________________________ read genotype  _____________________________________________

def read_genotype_data(file, sample_index, genotype_index, type_index, separator, header):
    '''Read and process the genotype file data, extracting the required columns into a pandas DataFrame.'''
    return read_data(file, [sample_index, genotype_index, type_index], ['Sample', 'Genotype', 'type'], separator, header)

# _____________________________________________ statistics on  the phenotype data  _____________________________________________
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



# QUANTITATIVE PHENOTYPE BINARISATION METHODS:
# _____________________________________________ #1 binariation  by median  _____________________________________________
@register_with_info("Average accuracy using median value for binarisation: ")
def binarisation_by_median(df, statistics, flip):
    '''Binarize the phenotype data based on the median.'''
    
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    if not flip:
        df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['median']).astype(int) # overwrite Phenotype values
    else:
        df_copy['Phenotype'] = (df_copy['Phenotype'] < statistics['median']).astype(int) # change the compare sign to obtain opposite values

    return df_copy

# _____________________________________________ #2 binariation  by extremes  _____________________________________________
# 2nd approach to binarisation
@register_with_info("Average accuracy using binarisation by extremes defined by 1st quartile and 3rd quartile: ")
def binarisation_by_extremes(df, statistics, flip):
    '''
    Binarize the phenotype data based on defined extreme intervals:
    - First extreme interval: values from the minimum to the first quartile are classified as 0.
    - Second extreme interval: values from the third quartile to the maximum are classified as 1.
    Values outside these intervals are set to -2.
    '''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original

    if not flip:
        choices = [0, 1] # Values to assign based on the conditions
    else:
        choices = [1, 0] # Values to assign based on the conditions

    conditions = [
    df_copy['Phenotype'] < statistics['first_quartile'], # zeros
    df_copy['Phenotype'] > statistics['third_quartile'] # ones
    ]

    # Apply conditions to the data, setting out-of-range values to -2
    df_copy['Phenotype'] = np.select(conditions, choices, default = -2)
    return df_copy

# _____________________________________________ #3 binariation  by 1st quartile  _____________________________________________
@register_with_info(f"Average accuracy using 1st quartile value for binarisation: ")
def binarisation_by_first_quartile(df, statistics, flip):
    '''Binarize the phenotype data based on the first quartile.'''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    
    if not flip:
        df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['first_quartile']).astype(int)
    else:
        df_copy['Phenotype'] = (df_copy['Phenotype'] < statistics['first_quartile']).astype(int)

    return df_copy

# _____________________________________________ #4 binariation  by 3rd quartile _____________________________________________
@register_with_info("Average accuracy using 3rd quartile value for binarisation: ")
def binarisation_by_third_quartile(df, statistics, flip):
    '''Binarize the phenotype data based on the third quartile.'''
    df_copy = df.copy() # Create a copy of the DataFrame to avoid modifying the original
    if not flip:
        df_copy['Phenotype'] = (df_copy['Phenotype'] > statistics['third_quartile']).astype(int)
    else:
        df_copy['Phenotype'] = (df_copy['Phenotype'] < statistics['third_quartile']).astype(int)

    return df_copy

# _____________________________________________ #5 all combinations for binariation  _____________________________________________
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

# _____________________________________________ average accuracy _____________________________________________

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

# _____________________________________________ histogram _____________________________________________
def histogram(df):
    ''' Create and plot the histogram for the phenotype data and analyze number of modes and the local maximums'''
    
    # We will also analyze the distribution, if the distribution is multimodal we want to detect modes and local maxims/ minims 
    def gaussian_kernel(u):
        '''
        Function to calculate Gaussian kernel for each element in array 'u'
        Gaussian Kernel formula: K(u) = (1 / sqrt(2π)) * exp(-0.5 * u^2)
        - 'u' represents the normalized distances between data points and the 
        evaluation point, scaled by the bandwidth.
        returns: An array of Gaussian kernel values corresponding to each distance in 'u'
        '''
        return (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * u**2)  # The Kernel function formula (Gaussian kernel is a common choice)
    
    def kernel_density_estimation(data, x_grid, bandwidth):
        '''
        Function to calculate the Kernel Density Estimation (KDE) over a set of points
        formula: f(x) = (1 / (n * h)) * Σ K((x - xi) / h)
        KDE formula: f(x) = (1 / (n * bandwidth)) * Σ K((x - xi) / bandwidth)
        returns: KDE values evaluated at each point in x_grid.
        '''
        u = (x_grid[:, None] - data) / bandwidth # 2D array of distances, each row represents a grid point and each column represents a data point
        contributions = gaussian_kernel(u) # 2D array, row represents a grid point, column represents a data point and each element represents an influence of data point to grid point
        kde_values = np.sum(contributions, axis=1) / (len(data) * bandwidth) # Sum these kernel values across data points (columns) for each grid point (rows).
        kde_values /= (np.sum(kde_values) * (x_grid[1] - x_grid[0])) # normalization to construct PDF (probability density function)
        return kde_values

    data = phenotype_df['Phenotype'].dropna().values
    
    # Example data generation
    # x = np.linspace(0, 10, 1000)
    # y = np.sin(x) * 100 + np.random.normal(0, 10, 1000)
    # data = y + np.random.normal(0, 30, 1000)

    # Preparation for KDE evaluations
    x_grid = np.linspace(data.min() - 1, data.max() + 1, 1000) # set of points for KDE evaluation (we cover the whole range of data)
    bandwidth= np.std(data) * np.power(4/(3*len(data)), 1/5)  # Silverman's rule of thumb for bandwidth selection; higher bandwidth = less sensitive, lower bandwidth = more sensitive to 'peaks' that do not represent modes but rather are minor fluctuations in data
    
    # KDE calculation
    kde_values = kernel_density_estimation(data, x_grid, bandwidth) # array of n elements

    # Histogram and KDE function calculations
    counts, bin_edges = np.histogram(data, bins='auto')  
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]
    kde_values = kde_values * len(data) * bin_width  # scale kde values, becaus it is a probability density function, but we represent histogram by counts of accesions in each bin not the probability

    # plot histogram and approximation line
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bin_edges, color='#93BFCF', edgecolor='#6096B4', histtype='stepfilled', linewidth=1.5, alpha=0.5, label='Histogram') # plotting histogram
    plt.plot(x_grid, kde_values, color='green', linewidth=3, alpha=0.3, label = 'scaled KDE to approximate histogram modes') # plotting the KDE function
    
    # Peak and modes detection using the first derivative of KDE
    kde_diff = np.diff(kde_values)  # First derivative of the KDE function, array of n-1 values representing differences
    # chceking where the sign of derivation changes (element-wise product of consecutive derivative values is negative and check if first value of pair was positive which means function was increasing before it started decreasing = peak )
    peaks = (kde_diff[:-1] * kde_diff[1:] < 0) & (kde_diff[:-1] > 0) # produces a boolean array where True symbolizes a detected peak
    peak_x = x_grid[1:-1][peaks]  # Align x_grid with peaks to obtain x indices of the peaks

    # Histogram peak refinement around KDE peaks
    search_radius = 2  # Adjusted search radius: number of bins to search around the KDE peak
    peak_indices = np.digitize(peak_x, bin_edges) - 1 # add the KDE x values for peaks to corresponding histogram bins
    refined_peaks = [] # array for the potential highest peaks in the modes
    for idx in peak_indices:
        local_indices = range(max(0, idx - search_radius), min(len(bin_centers), idx + search_radius + 1)) # bins around the detected peaks (the area ranges from detected peak -  search radius to detected peak +search radius)
        local_max = np.argmax(counts[local_indices]) # returns the index of highest bin in the searched area
        if counts[local_indices[local_max]] > max(counts)*0.03: # treshhold, dont consider peaks smaller than 3 % of the highest peak
            refined_peaks.append(local_indices[local_max]) # add the indices to the array that is later used for plotting the peak dot
        
    # Plot refined histogram peaks
    plt.scatter(bin_centers[refined_peaks], counts[refined_peaks], color='green', label='Histogram Peaks')

    # Add grid, labels, and title
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Phenotype Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of trait values')
    plt.savefig("histogram.png", dpi=300, bbox_inches='tight')
    plt.show()
    print(f'Detected number of modes in histogram is: {len(refined_peaks)} with the peak frequency values:{counts[refined_peaks]}.')

# _____________________________________________ boxplot _____________________________________________
def boxplot(genotype_df, phenotype_df):
    '''
    Visualize the distribution of phenotype values for accessions with 
    reference (0) and alternative (1) genotypes.
    '''
    merged_df = pd.merge(genotype_df, phenotype_df, on='Sample')

    # Create a figure and a set of subplots
    fig, ax = plt.subplots(figsize=(15, 10))

    # Predefining required properties for boxplot: the box and median line, (whiskers and caps can be added too)
    box_properties = dict(linestyle='-', linewidth=2, color='#6096B4', facecolor='#93BFCF', alpha=0.5)
    median_properties = dict(linestyle='-', linewidth=2.5, color='firebrick')
    # whisker_properties = dict(linestyle='-', linewidth=3, color='#6097B466')
    # caps_properties = dict(linestyle='-', linewidth=3, color='#6097B466')

    # Create the boxplot with the properties for box and median
    # 'patch_artist' must be set on True for properties like facecolor and alpha to work
    bp = merged_df.boxplot(
        column='Phenotype',
        by='Genotype',
        ax=ax,
        boxprops=box_properties, 
        medianprops=median_properties,
        showfliers=False,
        showcaps=False,
        whis = 0,
        patch_artist=True,
        )
    
    # use another color for each type of accesions
    conditions = [
                (merged_df['type'] == 'Elite'),
                (merged_df['type'] == 'G. soja'),
                (merged_df['type'] == 'Genetic'),
                (merged_df['type'] == 'Landrace')
                ]
    choices = ['blue', 'green', 'red', 'darkviolet']  # Corresponding values for the types above
    
    merged_df['color'] = np.select(conditions, choices, default = 'black')
    
    # Adding jitter and  generating correct boxplots columns labels
    # create 2 groups based on genotype 
    # - Group 1: reference genotype rows
    # - Group 2: alternative genotype rows
    # enumerate encapsulates all objects (DataFrameGroupBy) merged_df.groupby returns (In this case 2 objects: Both objects contain two values (key, group)
    name_setter = [] # list for labels describing a boxplot ('reference', 'alternative')
    name_setter_index = [] # correct indexes for titles, boxplots are indexed from 1, because index 0 belongs to y ax
    counts_per_group = []  # Store counts to display later
    for i, (key, group) in enumerate(merged_df.groupby('Genotype')):
         # Append labels and indices
        if key == 0:
            name_setter.append(f'REF Genotype ({len(group)} accesions)')
            name_setter_index.append(i + 1)
        elif key == 1:
            name_setter.append(f'ALT Genotype ({len(group)}accesions)')
            name_setter_index.append(i + 1)
        counts_per_group.append(len(group))

        # jitered x values:
        x = np.random.normal(i + 1, 0.02, size=len(group))
        ax.scatter(x, group['Phenotype'], alpha=0.5, color=group['color'], edgecolor='none')

    #adding a legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Elite', markerfacecolor='blue', alpha = 0.5, markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='G. soja', markerfacecolor='green', alpha = 0.5, markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Genetic', markerfacecolor='red', alpha = 0.5, markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Landrace', markerfacecolor='purple', alpha = 0.5, markersize=10)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize='large')

    # additional plot formatting
    ax.set_title('Phenotype Distribution for accesions with reference and alternative genotype')
    plt.suptitle('')  # Removes the 'Group By' title
    plt.xlabel('')
    plt.ylabel('Phenotype Value')
    plt.xticks(name_setter_index, name_setter)  # Sets custom x-axis labels

    # Show the plot
    plt.savefig("boxplot.png", dpi=300, bbox_inches='tight')
    plt.show()

# _____________________________________________ MAIN _____________________________________________

# CONSTANT VARIABLES; adjust for your specific phenotype / genotype files characteristics
phenotype_file = 'aggregated_mean.csv'
phenotype_file_header_present = 'Y'
phenotype_file_sampleID_column_index = 0 # index of accesion numbers column in phenotype file (indexing starts from 0)
phenotype_file_values_column_index = 1 # index of phenotype values column in phenotype file (indexing starts from 0)

genotype_file = 'genotype_binarised_data.csv'
genotype_file_header_present = 'Y'
genotype_file_sampleID__column_index = 1 # index of accesion numbers column in genotype file (indexing starts from 0)
genotype_file_values_column_index = 2 # index of binarised genotype values column in genotype file (indexing starts from 0)
genotype_file_type_column_index = 3 # index of type column in genotype file (landrace, cultivar ...)

separator = ',' # both genotype and phenotype file should use same separator
flip = True # change to True if needed


# load the data from phenotype and genotype file into dataframe

phenotype_df = read_phenotype_data(phenotype_file, phenotype_file_sampleID_column_index, phenotype_file_values_column_index, separator, phenotype_file_header_present)
binarised_genotype_df = read_genotype_data(genotype_file, genotype_file_sampleID__column_index, genotype_file_values_column_index, genotype_file_type_column_index, separator, genotype_file_header_present)


stats = compute_statistics(phenotype_df)

for function, info in zip(function_list, info_strings):
    binarised_phenotype_df = function(phenotype_df, stats, flip)
    avg = average_accuracy_from_dataframes(binarised_genotype_df, binarised_phenotype_df)
    print(f'{info}{avg:.2f}%')

binarisation_combinations(binarised_genotype_df, phenotype_df)

# checking how many genotype file accesions do not have match in phenotype file
mask = ~binarised_genotype_df['Sample'].isin(phenotype_df['Sample'])
result_df = binarised_genotype_df[mask]
accesions_without_phenotype = len(result_df)

# OUTPUT
histogram(phenotype_df)
boxplot(binarised_genotype_df, phenotype_df)

print(f'Number of accesions that could not be matched to phenotype value: {accesions_without_phenotype}.')
print(f"The phenotype file contains {len(phenotype_df.Phenotype)} samples.")
print(f"The MINIMAL VALUE of the trait is:   {stats['min']:.2f}.")
print(f"The MAXIMAL VALUE of the trait is:   {stats['max']:.2f}.")
print(f"The AVERAGE VALUE of the trait is:   {stats['mean']:.2f}.")
print(f"The MEDIAN VALUE of the trait is:    {stats['median']:.2f}.")
print(f"The 1ST QUARTILE of the trait is:    {stats['first_quartile']:.2f}.")
print(f"The 3RD QUARTILE of the trait is:    {stats['third_quartile']:.2f}.")