import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], delimiter='\t')
df['GC content'] = df.apply(lambda row: (row['RefSeq'].count('G')+row['RefSeq'].count('C'))/len(row['RefSeq']), axis=1) 
tmp = df.groupby(['RefSeq', 'series']).size()

tmp.to_csv(sys.argv[2])
