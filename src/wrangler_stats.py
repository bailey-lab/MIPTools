import pandas as pd

def strip_percent(value):
    try:
        value=int(str(value).split("(")[0])
    except (ValueError, TypeError):
        pass
    return value

def get_stats(stat_file):
    table=pd.read_csv(stat_file, sep='\t')
    header=list(table.columns.values)
    for column in header:
        table[column]=table[column].map(strip_percent)
    return table