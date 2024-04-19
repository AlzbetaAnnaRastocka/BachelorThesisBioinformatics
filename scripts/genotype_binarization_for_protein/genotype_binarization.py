import csv

# output file
file = 'Soy2939_Glyma.20g085100_All_Data.csv'
output_file_name = 'genotype_binarised_data.csv'

#open input file that we want to filter
with open(file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    # read the header row of the input file to get the column names
    headers = next(csv_reader)
    # create a dictionary to store indexes for columns in the header
    header_map = {name: index for index, name in enumerate(headers)}

    # open output file for the filtered data
    with open('genotype_binarised_data.csv', 'w', newline='') as output_file:
        csv_writer = csv.writer(output_file, delimiter=',')
        
        # header of output file
        csv_writer.writerow(['Accesion', 'GRIN_Accesion', 'WT/MUT', 'type'])
        
        # loop through the rows in the input file
        for row in csv_reader:
            # check for the target string (mutant) in the "Genotype_Description" column
            if "conservative_inframe_insertion|" in row[header_map['Genotype_Description']]:  # column index 12 is the 13th column
                # write only the "Accesion" column and "GRIN_Accesion" column along with the value 1 (mutant) in the WT/MUT column.
                csv_writer.writerow([row[header_map['Accession']], row[header_map['GRIN_Accession']], 1, row[header_map['Improvement_Status']]])

            # check for the target string (reference genome) in the "Genotype_Description" column
            else:
                csv_writer.writerow([row[header_map['Accession']], row[header_map['GRIN_Accession']], 0, row[header_map['Improvement_Status']]]) #value 0 (wild type)
