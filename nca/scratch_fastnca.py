### create NCA analysis
import json
import pandas as pd 
#df=pd.read_excel('/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/Supplementary_Table_1_Trancriptional_reporter_library_screening.xlsx')
# Load the specific sheets from the Excel file
file_path = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/Supplementary_Table_1_Trancriptional_reporter_library_screening.xlsx'

# Read each sheet into separate DataFrames
df_screening_lb = pd.read_excel(file_path, sheet_name='Screening_LB')
df_screening_m9 = pd.read_excel(file_path, sheet_name='Screening_M9')

# Display the head of each DataFrame to confirm
print("Screening LB DataFrame:")
print(df_screening_lb.head())

print("\nScreening M9 DataFrame:")
print(df_screening_m9.head())

#### get crp regulation kegg enrichment analysis
# Replace 'file.txt' with your file path
crp_m9 = pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/m9_crp.txt', header=None)
crp_lb = pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/lb_crp.txt', header=None)





####

json_file='/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/regulatoryNetwork_RDBECOLITFC00189.json'

# Parse the JSON data
# Read and parse JSON data from a file
with open(json_file, 'r') as file:
    parsed_data = json.load(file)
# Extract node information
bolArray=[]
nodes = []
genes=[]
for node in parsed_data["elements"]["nodes"]:
    node_info = {
        "id": node["data"]["id"],
        "label": node["data"]["label"],
        "network": node["data"]["network"],
        "isChild": node["data"]["isChild"],
        "isParent": node["data"]["isParent"],
        "effect": node["data"]["effect"],
        "position": node["position"],
        "group": node["group"]
    }
    nodes.append(node_info)
    if node_info['isChild'] and (node_info['effect']=='activator'):
        genes.append(node_info['label'])
        bolArray.append(1)
    elif node_info['isChild'] and (node_info['effect']=='repressor'):
        genes.append(node_info['label'])
        bolArray.append(-1)
    else:
        print(node_info)



       
df=pd.DataFrame()
df['genes']=genes
df['CRP']=bolArray

# Filter df_screening_lb to only include genes in df['genes']
filtered_screening = df_screening_lb[df_screening_lb['Gene_Name'].isin(df['genes'])]

# Merge to get common genes in the same order as in df['genes']
result = pd.merge(df, filtered_screening, left_on='genes', right_on='Gene_Name', how='left')
# Replace NaN values with 0
result = result.fillna(0)


# Filter df_screening_lb to only include genes in df['genes']
filtered_screening = df_screening_m9[df_screening_m9['Gene_Name'].isin(df['genes'])]

# Merge to get common genes in the same order as in df['genes']
result2 = pd.merge(df, filtered_screening, left_on='genes', right_on='Gene_Name', how='left')

# Replace NaN values with 0
result2 = result2.fillna(0)

### get columns
cols=['genes','CRP','ratio_t0', 'ratio_t5', 'ratio_t30',
       'ratio_t60', 'ratio_t120', 'ratio_t180']

df_m9_crp=result.copy()
df_lb_crp=result2.copy()


df_m9_crp.to_excel('crp_net_data_m9.xlsx')
df_lb_crp.to_excel('crp_net_data_lb.xlsx')

import pdb;pdb.set_trace()
# Print extracted nodes
# for node in nodes:
#     print(node)


