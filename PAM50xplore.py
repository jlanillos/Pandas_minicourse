# Dataframes: https://drive.google.com/drive/folders/1luAZrsdESz_r6Xi3iy9eyhiTxzg5x9LM?usp=sharing
import pandas as pd
import numpy as np


# rnacounts is the dataframe with the counts from TCGA-BRCA rnaseq per gene (Xena Browser)
rnacounts = pd.read_csv('TCGA-BRCA.htseq_counts.tsv',sep='\t') # Open rnaseq counts matrix and load it into a dataframe called rnacount


# Here, let's assume this is the list of 50 genes related to PAM50 calculations
pam50 = 'ACTR3B,ANLN,BAG1,BCL2,BIRC5,BLVRA,CCNB1,CCNE1,CDC20,CDC6,CDH3,CENPF,CEP55,CXXC5,EGFR,ERBB2,ESR1,EXO1,FGFR4,FOXA1,FOXC1,GPR160,GRB7,KIF2C,KRT14,KRT17,KRT5,MAPT,MDM2,MELK,MIA,MKI67,MLPH,MMP11,MYBL2,MYC,NAT1,NDC80,NUF2,ORC6,PGR,PHGDH,PTTG1,RRM2,SFRP1,SLC39A6,TMEM45B,TYMS,UBE2C,UBE2T'
pam50 = pam50.split(',') # Convert the string variable pam50 to a list splitting by commas (',')
len(pam50) # Let's check the length of the list "pam50"


# Let's get the dataframe with Ensembl to gene symbol convertion
appris = pd.read_csv('appris_data.appris_label.csv',sep='\t') # Open a table with the ensembl gene IDs to gene synbol convertion
appris = appris[['ensembl_id','SYMBOL']] # We only need the data from the first two columns --> ['ensembl_id',''SYMBOL']
appris.head() # Let's print the first rows of this dataframe using .head() function
appris = appris.drop_duplicates(subset='SYMBOL', keep='first') # Let's drop duplicated line using .drop_duplicates() function
appris.head() # Let's print the first rows of appris dataframe using .head() function
appris = appris[appris['SYMBOL'].isin(pam50)] # Let's filter only the 50 genes of interest using .isin() function
len(appris) # Let's check the resulting length (number of rows) in the appris dataframe after filtering only for the 50 genes of interest
appris.head() # Check the name of the columns we really need in the appris table


# Learning checkpoint: load a table into a dataframe (also .read_excel(), display few rows of a dataframe, get specific columns from a dataframe

ensembl_genes = list(appris['ensembl_id']) # create a list with the gene symbols in the first column of appris dataframe


# Let's find the 50 genes of interest in the rnacounts dataframe and discard the rest of rows
rnacounts.head() # Display first rows of the dataframe to recall column names we need
rnacounts['Ensembl_ID'] = rnacounts['Ensembl_ID'].str.split('.').str[0] # Remove ensembl_ID version
rnacounts.head() # Display first rows of the dataframe to ckech the change in "Ensembl_ID" column
len(rnacounts)
rnacounts = rnacounts.loc[rnacounts['Ensembl_ID'].isin(ensembl_genes)] # usethe .isin function to find rows matching genes of interest in the "Ensembl_ID" column
len(rnacounts)


# Learning checkpoint: apply a "string" operation over a column, select rows searching multiple values in a specific column


# Change the "ensembl" ids by gene "symbol" ids in the rnacounts matrix
ensembl_genes = list(appris['ensembl_id']) # the keys of the dictionary # We have already execute this line above, it's a reminder
gene_symbols = list(appris['SYMBOL']) # the values of the dictionary
d = zip(ensembl_genes,gene_symbols) # Create the dictionary from appris columns (part 1: zip two lists [keys,values])
d = dict(d) # Create the dictionary from appris columns (part 2: Convert the "zipped" object into the dictionary)
print(d)
rnacounts['SYMBOL'] = rnacounts['Ensembl_ID'].map(d) # Using .map() function, we translate "ensembl" IDs to "symbol" IDs (this is as useful as the vlookup table in Excel/Libreoffice)
rnacounts.head()
rnacounts = rnacounts.drop(columns=['Ensembl_ID']) # Let's delete "Ensembl_ID' column (using .drop() function) 
rnacounts.head()
rnacounts = rnacounts.set_index('SYMBOL') # Let's set "SYMBOL" column as the dataframe index (using .set_index() function)
rnacounts.head()


# Learning checkpoint: dictionary from two lists and .map() the dictionary (similarly as vlookup in Excel) to a new column, delete colums (.drop())



##### Question: Does a PCA analysis (2D) of TCGA-BRCA RNAseq counts recapitulate PAM50 classification? #####


# Let's visualize a PCA (2D) plot of the PAM50 genes using the rnacounts
# Learning goals: import multiple Python libraries and transform the dataframe: transpose and apply a fitting
# Learning goals: introduce the .loc() function
# Learning goals: practice some learning outcomes from above (load dataframes, map dictionaries)
# Learning goals: scatter plot using Pandas and matplotlib

# Load some libraries and utilities
from sklearn.preprocessing import StandardScaler # we need this library to "normalize/transform" our data
from sklearn.decomposition import PCA # we need this library to perform PCA
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


x = rnacounts.T # Let's transpose the rnacounts table and save it into x
x.head()
samples = list(x.index)

x = StandardScaler().fit_transform(x) # we are fitting (applying a transformation) in our data
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)

PCA_df = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
PCA_df['samples'] = samples
PCA_df.head()


# Let's add PAM50 classification 
pheno = pd.read_csv('brca.pam50.csv',sep='\t') # Load a dataframe with PAM50 classification of a number of samples, how many?
pheno.head()
len(pheno)
pam50samples = list(pheno['sample'])
pam50vals = list(pheno['PAM50'])
d = dict(zip(pam50samples,pam50vals)) # Again, create a dictionary with the keys and values that we will map in the matrix
print(d)
PCA_df['PAM50'] = PCA_df['samples'].map(d)
# Setting samples as row "indexes" of the dataframe
PCA_df = PCA_df.set_index('samples')
PCA_df.head()
# Assign 'Not reported value in "PAM50" column to those samples with no PAM50 classification
PCA_df.loc[(PCA_df['PAM50'].isnull()),'PAM50']='Not reported'


# Create dataframes df2plot and df2plot_reported
df2plot = PCA_df
colors = ['purple','orange','cyan','red', 'green', 'blue','gray'] # A list with as many colors as the number of "PAM50" subtypes
targets = list(set(df2plot['PAM50'])) # ['Unknown', 'Basal', 'Normal', 'LumA', 'Her2', 'LumB', 'Not reported'] #PAM50" subtypes"
color_dict = dict(zip(targets,colors)) # Again, create a dictionary to map the correct color to each sample based on their "PAM50" classification
df2plot['color'] = df2plot['PAM50'].map(color_dict)

df2plot_reported = df2plot.loc[df2plot['PAM50'] != 'Not reported'] # Let's filter out samples with no (NaN) PAM50 label using .loc and output the result into a new dataframe 



# PLot PCA (only samples with PAM50 info)
fig, ax1 = plt.subplots(1,1,figsize=(10,8)) # Creating a "Figure object"
df2plot_reported.plot.scatter(ax=ax1,x='principal component 1', y='principal component 2', c='color', s = 20, colorbar=False) # Plotting a scatter from a pandas dataframe

# Add some info to the figure
ax1.set_xlabel('Principal Component 1', fontsize = 12)
ax1.set_ylabel('Principal Component 2', fontsize = 12)
fig.suptitle('TCGA-BRCA RNAseq counts 2 component PCA', fontsize = 14)
ax1.set_title('Samples with PAM50 classification')
patches = list()
for category in color_dict.keys():
    patches.append(mpatches.Patch(color=color_dict[category], label=category))
ax1.legend(handles=patches, bbox_to_anchor=(1.01, 0.5), loc='lower left')
#840
ax1.grid()

#Display the plot
plt.show()



# We are going to repeat the plot, but instead just displaying the samples with PAM50 classification, let's show all samples in the TCGA-BRCA cohort, and see how they are distributed

# PLot PCA (also not reported)
fig,(ax1, ax2) = plt.subplots(1,2,figsize=(17,7))

# There will be two plots. The one on the left, as previously done, it will only display PAM50-annotated samples
df2plot_reported.plot.scatter(ax=ax1,x='principal component 1', y='principal component 2', c='color', s = 20, colorbar=False)

# The second plot will be on the right, and will show all TCGA-BRCA samples
df2plot.plot.scatter(ax=ax2,x='principal component 1', y='principal component 2', c='color', s = 20, alpha = 0.5, label='PAM50', colorbar=False)


ax1.set_xlabel('Principal Component 1', fontsize = 15)
ax2.set_xlabel('Principal Component 1', fontsize = 15)
ax1.set_ylabel('Principal Component 2', fontsize = 15)
fig.suptitle('TCGA-BRCA RNAseq counts 2 component PCA', fontsize = 20)
ax1.set_title('Samples with PAM50 classification')
ax2.set_title('All samples')
patches = list()
for category in color_dict.keys():
    patches.append(mpatches.Patch(color=color_dict[category], label=category))
ax2.legend(handles=patches, bbox_to_anchor=(1.01, 0.5), loc='lower left')
ax1.grid()
ax2.grid()

plt.show()




#### How to save dataframes into CSV files or Excel files?

# In .csv format
df2plot.to_csv('PCA_analysis_TCGA-BRCA_PAM50.csv',sep='\t')
# In .xlsx format
df2plot.to_excel('PCA_analysis_TCGA-BRCA_PAM50.xlsx')





################################################################
########################### BONUS ##############################
################################################################

### Bonus: .groupby()
pheno.head()
pheno.groupby('PAM50').count()
pheno.groupby('PAM50').count().to_csv('PAM50_summary.csv',sep='\t')


## Bonus: stats
rnacounts.head()

# Apply basic operations in one sample
rnacounts['TCGA-E9-A1NI-01A'].mean()
rnacounts['TCGA-A7-A13F-01A'].std()
rnacounts['TCGA-E9-A1NI-01A'].sum()
# Or all samples (the whole dataframe)
rnacounts.mean()

## Stats to copare sets of samples:
pheno.head()
pheno['sample'].loc[pheno['PAM50'] == 'LumA']

lumA_samples =  list(pheno['sample'].loc[pheno['PAM50'] == 'LumA'])
Basal_samples =  list(pheno['sample'].loc[pheno['PAM50'] == 'Basal'])

# Combination of the .apply() function + "lambda" variable (called x in this example) to make a t-test (using scipy.stats.ttest_ind function) between LuminalA and Basal samples
rnacounts['ttest_Luminal'] = rnacounts.apply(lambda x: scipy.stats.ttest_ind(x[lumA_samples], x[Basal_samples]), axis = 1)



