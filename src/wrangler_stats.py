
import pandas as pd


def remove_percent(item):
    try:
        return int(item.split("(")[0])
    except (IndexError, ValueError):
        return item


def get_stats(stat_file):
    sti = pd.read_table(stat_file).applymap(remove_percent)
    sti = sti.loc[sti["Sample"] != "Sample"]
    for c in sti.columns:
        try:
            sti[c] = sti[c].astype(int)
        except ValueError:
            pass
    return sti
