import pandas as pd
import pingouin as pg

## Load data
datafile = "tmp.csv" 
df = pd.read_csv(datafile)

## Compute ANOVA
paovm=pg.rm_anova(data=df,dv='fraction_correct',within=['presentation_condition','source_condition'],subject='subject',correction='auto',detailed=True,export_filename='tmpStts')
print("=== Pingouin ANOVA === sphericity: ", pg.sphericity(paovm))
pg.print_table(paovm)

