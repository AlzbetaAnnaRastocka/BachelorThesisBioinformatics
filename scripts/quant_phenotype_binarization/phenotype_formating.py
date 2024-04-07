import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Original phenotype file name
file = "aggregated_mean.csv"
# Original phenotype file phenotype column number (n-1, for phenotype in 2nd column use number 1)
phenotype_column_index = 1
# Original phenotype file Accession number column number (n-1, for phenotype in 2nd column use number 1)
sample_ID_index = 0
# Original phenotype file  column separator (examples: "\t",",",";")
separator = ","
# Original phenotype file sample name format (G - Grin name format, S - Soy1066 name format , S7 -Soy775 name format )
namform = "G"
# Original phenotype type (Qual - qualitative, Quant - quantitative)
phentyp = "Quant"
#Header in the Original phenotype file ("Y","N")
header="Y"

#Result phenotype file name
res = "pheno_res.txt"
#Result phenotype file formate ("AccuCalc","MADis","AccuTool")
form = "AccuCalc"


# creating two dataframes (Phenotype and Sample)

file = open(file,"rt")

name_list,phen_list = list(),list()
n = 0

for line in file:
    if header == "Y" and n == 0: n = 1
    else:
        line=line.split(separator)
        name_list.append(line[sample_ID_index])
        phen_list.append(line[phenotype_column_index])

file.close()

df = pd.DataFrame()
df["Sample"] = name_list
df["Phenotype"] = phen_list


# Data statistic
print (f"The phenotype file contains {len(df.Phenotype)} samples.")

if phentyp == "Quant":
    df["Phenotype"]=pd.to_numeric(df["Phenotype"],errors='coerce')
    first_quartile = df.Phenotype.quantile(0.25)
    second_quartile = df.Phenotype.quantile(0.5)
    third_quartile = df.Phenotype.quantile(0.75)
    max_value = df.Phenotype.max()
    min_value = df.Phenotype.min()
    mean_value = df.Phenotype.mean()
    delt = max_value - min_value

    # # phenotype binarisaiton

    # 1st Approach
    df['WT/MUT'] = (df['Phenotype'] > second_quartile).astype(int)

    df[['Sample', 'WT/MUT']].to_csv('median.csv', index=False)

    # 2nd Approach
    conditions = [
        df['Phenotype'] < first_quartile, # zeros - we take 0-25% of values
        df['Phenotype'] > third_quartile # ones - we take 75-100% of values
    ]
    choices = [0, 1]  # Corresponding values for the conditions above

    # we take conditions and apply choices to "Extremes" by those conditions, 
    # if no record in condition, we set value to NaN
    df['WT/MUT'] = np.select(conditions, choices, default = -2) 

    #df['WT/MUT'] = df['WT/MUT'].astype(int) # convert float type to int
    df[['Sample', 'WT/MUT']].to_csv('extremes.csv', index=False) 

    # 3rd Approach
    df['WT/MUT'] = (df['Phenotype'] > first_quartile).astype(int)
    df[['Sample', 'WT/MUT']].to_csv('first_quartile.csv', index=False)
    # 4th Approach
    df['WT/MUT'] = (df['Phenotype'] > third_quartile).astype(int)
    df[['Sample', 'WT/MUT']].to_csv('third_quartile.csv', index=False)

    print (f"The minimal value of the trait is {min_value} and the maximal value is {max_value}.")
    print (f"The average of the trait is {mean_value} and the median value is {second_quartile}.")
    print (f"The 25% quantile (1st) of the trait is {first_quartile} and the 75% quantile (3rd) of is {third_quartile}.")

    if delt <3:plt.hist(df.Phenotype, bins = 20)
    elif delt <5:plt.hist(df.Phenotype, bins = 30)
    elif delt <7:plt.hist(df.Phenotype, bins = 50)
    elif delt <9:plt.hist(df.Phenotype, bins = 70)
    else:plt.hist(df.Phenotype, bins = 100)
    plt.savefig("histogram.png", dpi=300, bbox_inches='tight')
    plt.show()
else: 
    raise Exception("Incorrect phenotype format type.")


# average accuracy 
def average_accuracy(genotype_file, phenotype_file, separator=',', header=True):
    if header: header_val = 0 
    else: header_val = None
    
    # load csv files (phenotype, genotype) into dataframe
    df_genotype = pd.read_csv(genotype_file, sep=separator, header=header_val, names=["Soy2939", "GRIN_Accession", "Genotype"])
    df_phenotype = pd.read_csv(phenotype_file, sep=separator, header=header_val, names=["GRIN_Accession", "Phenotype"])
    
    # Merge the dataframes on the 'GRIN_Accession' column
    merged_df = pd.merge(df_genotype, df_phenotype, on="GRIN_Accession")
    merged_df.to_csv('merged.csv', index=False)

    # formulas needed for average accuracy
    N_WTcorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 0)])
    N_WTincorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["Phenotype"] == 1)])
    acc_WT = (N_WTcorrect / (N_WTcorrect + N_WTincorrect)) * 100

    N_MUTcorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 1)])
    N_MUTincorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["Phenotype"] == 0)])
    acc_MUT = (N_MUTcorrect / (N_MUTcorrect + N_MUTincorrect)) * 100

    # average accuracy formula
    return (acc_MUT + acc_WT) / 2
    
avg = average_accuracy("filtered_data.csv", "median.csv")
print (f"Average accuracy using median method for phenotype binarisation is: {avg}.")

avg = average_accuracy("filtered_data.csv", "extremes.csv")
print (f"Average accuracy using extremes method for phenotype binarisation is: {avg}.")

avg = average_accuracy("filtered_data.csv", "first_quartile.csv")
print (f"Average accuracy using first_quartile method for phenotype binarisation is: {avg}.")

avg = average_accuracy("filtered_data.csv", "third_quartile.csv")
print (f"Average accuracy using third_quartile method for phenotype binarisation is: {avg}.")