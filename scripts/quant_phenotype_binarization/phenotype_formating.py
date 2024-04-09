import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Original phenotype phenotype_file name
phenotype_file = "aggregated_mean.csv"
genotype_file = "filtered_data.csv"

# Original phenotype phenotype_file phenotype column number (n-1, for phenotype in 2nd column use number 1)
phenotype_column_index = 1
# Original phenotype phenotype_file Accession number column number (n-1, for phenotype in 2nd column use number 1)
sample_ID_index = 0
# Original phenotype phenotype_file  column separator (examples: "\t",",",";")
separator = ","
# Original phenotype phenotype_file sample name format (G - Grin name format, S - Soy1066 name format , S7 -Soy775 name format )
namform = "G"
# Original phenotype type (Qual - qualitative, Quant - quantitative)
phentyp = "Quant"
#Header in the Original phenotype phenotype_file ("Y","N")
header="Y"

#Result phenotype phenotype_file name
res = "pheno_res.txt"
#Result phenotype phenotype_file formate ("AccuCalc","MADis","AccuTool")
form = "AccuCalc"


# creating two dataframes (Phenotype and Sample)

phenotype_file = open(phenotype_file,"rt")

name_list,phen_list = list(),list()
n = 0

for line in phenotype_file:
    if header == "Y" and n == 0: n = 1
    else:
        line=line.split(separator)
        name_list.append(line[sample_ID_index])
        phen_list.append(line[phenotype_column_index])

phenotype_file.close()

df = pd.DataFrame()
df["Sample"] = name_list
df["Phenotype"] = phen_list


# Data statistic
print (f"The phenotype phenotype_file contains {len(df.Phenotype)} samples.")

if phentyp == "Quant":
    df["Phenotype"]=pd.to_numeric(df["Phenotype"],errors='coerce')
    first_quartile = df.Phenotype.quantile(0.25)
    second_quartile = df.Phenotype.quantile(0.5)
    third_quartile = df.Phenotype.quantile(0.75)
    max_value = df.Phenotype.max()
    min_value = df.Phenotype.min()
    mean_value = df.Phenotype.mean()
    delt = max_value - min_value

    twenty = df.Phenotype.quantile(0.01)
    thirty = df.Phenotype.quantile(0.40)
    seventy = df.Phenotype.quantile(0.60)
    eighty = df.Phenotype.quantile(0.99)

    # phenotype binarisaiton

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

    # 5th approach
    conditions = [
        (df['Phenotype'] > twenty) & (df['Phenotype'] < thirty),  # zeros - we take 20 - 30% of values
        (df['Phenotype'] > seventy) & (df['Phenotype'] < eighty)  # ones - we take 75-100% of values
    ]

    choices = [0, 1]  # Corresponding values for the conditions above

    df['WT/MUT'] = np.select(conditions, choices, default= "NaN")
    df[['Sample', 'WT/MUT']].to_csv('skuska.csv', index=False)




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
    
avg = average_accuracy(genotype_file, "median.csv")
print (f"Average accuracy using median method for phenotype binarisation is: {avg}.")

avg = average_accuracy(genotype_file, "extremes.csv")
print (f"Average accuracy using extremes method for phenotype binarisation is: {avg}.")

avg = average_accuracy(genotype_file, "first_quartile.csv")
print (f"Average accuracy using first_quartile method for phenotype binarisation is: {avg}.")

avg = average_accuracy(genotype_file, "third_quartile.csv")
print (f"Average accuracy using third_quartile method for phenotype binarisation is: {avg}.")

avg = average_accuracy(genotype_file, "skuska.csv")
print (f"Average accuracy using skuska method for phenotype binarisation is: {avg}.")


phenotype_file = open(genotype_file,"rt")

n_list, g_list = list(), list()
for line in phenotype_file:
        line=line.split(separator)
        line[2].strip()
        if line[2] == "0\n":
            n_list.append(line[1])
            g_list.append(line[2].strip())

phenotype_file.close()

dtf = pd.DataFrame()
dtf["GRIN_Accession"] = n_list
dtf["refgef"] = g_list
dtf[['GRIN_Accession', 'refgef']].to_csv('ref.csv', index=False)


dtf_phenotype = pd.read_csv("aggregated_mean.csv", sep = separator, header = 0, names = ["GRIN_Accession", "Phenotype"])

merged_dtf = pd.merge(dtf[['GRIN_Accession', 'refgef']], dtf_phenotype, on="GRIN_Accession")
merged_dtf.to_csv('joooj.csv', index=False)

max_value = merged_dtf.Phenotype.max()
min_value = merged_dtf.Phenotype.min()
print(max_value)
print(min_value)


#       skuska
def average_accuracy_from_dataframe(genotype_df, phenotype_df):
    # Merge the dataframes on the 'GRIN_Accession' column
    merged_df = pd.merge(genotype_df, phenotype_df, on="GRIN_Accession")

    # formulas needed for average accuracy
    N_WTcorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["WT/MUT"] == 0)])
    N_WTincorrect = len(merged_df[(merged_df["Genotype"] == 0) & (merged_df["WT/MUT"] == 1)])
    acc_WT = (N_WTcorrect / (N_WTcorrect + N_WTincorrect)) * 100

    N_MUTcorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["WT/MUT"] == 1)])
    N_MUTincorrect = len(merged_df[(merged_df["Genotype"] == 1) & (merged_df["WT/MUT"] == 0)])
    acc_MUT = (N_MUTcorrect / (N_MUTcorrect + N_MUTincorrect)) * 100

    # average accuracy formula
    return (acc_MUT + acc_WT) / 2

def binarisation_combinations (phenotype_csv):
    phenotype_dataframe = pd.read_csv(phenotype_csv, sep = separator, header = 0, names = ["GRIN_Accession", "Phenotype"])
    combinations_for_binarisation_dict = dict()
    interval = 50
    dictionary_index = 0
    for segment_zeros in range(0, 100, interval):
        for segment_ones in range(0, 100, interval):
            if segment_zeros == segment_ones:
                continue
            else:
                zeros_quantil = segment_zeros / 100
                zeros_quantil_end = (segment_zeros + interval) / 100
                ones_quantil = segment_ones / 100
                ones_quantil_end = (segment_ones + interval) / 100

                conditions = [
                    (phenotype_dataframe['Phenotype'] > phenotype_dataframe.Phenotype.quantile(zeros_quantil)) & (phenotype_dataframe['Phenotype'] < phenotype_dataframe.Phenotype.quantile(zeros_quantil_end)),
                    (phenotype_dataframe['Phenotype'] >  phenotype_dataframe.Phenotype.quantile(ones_quantil)) & (phenotype_dataframe['Phenotype'] < phenotype_dataframe.Phenotype.quantile(ones_quantil_end)) # ones - we take 75-100% of values
                    ]
                choices = [0, 1]  # Corresponding values for the conditions above

                phenotype_dataframe['WT/MUT'] = np.select(conditions, choices, default = -2) # check if correct data are saved

                dataframe_phenotype_copy = phenotype_dataframe[["GRIN_Accession", "WT/MUT"]].copy()
                # if dictionary_index % 10 == 0:
                #         dataframe_phenotype_copy.to_csv(f'ukazka{dictionary_index}.csv', index=False)


                # Forming the string
                info_string = f"Borders of 0s: {zeros_quantil}, {zeros_quantil_end} \nborders of 1s: {ones_quantil}, {ones_quantil_end}"

                # Updating the dictionary with a tuple containing the dataframe slice and the string
                combinations_for_binarisation_dict.update({ dictionary_index: (dataframe_phenotype_copy, info_string)})
                dictionary_index += 1

                #print(f"Borders of 0s: {zeros_quantil}, {zeros_quantil_end} \nborders of 1s: {ones_quantil}, {ones_quantil_end}\n")
    genotype_dataframe = pd.read_csv(genotype_file, sep=separator, header = 0, names=["Soy2939", "GRIN_Accession", "Genotype"])

    accuracies = dict()
    for key, value in combinations_for_binarisation_dict.items():
        # value is a tuple where value[0] is the dataframe and value[1] is the string
        accuracy = average_accuracy_from_dataframe(genotype_dataframe, value[0])
        accuracies.update({key: accuracy})

    key_with_highest_value = max(accuracies, key=accuracies.get)
    value_best = combinations_for_binarisation_dict[key_with_highest_value]
    print(f'Key with the highest value: {key_with_highest_value}, Value: {accuracies[key_with_highest_value]}, {value_best[1]}')


binarisation_combinations("aggregated_mean.csv")