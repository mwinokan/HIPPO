
import mout
import numpy as np

def df_row_to_dict(df_row):

    assert len(df_row) == 1

    data = {}

    for col in df_row.columns:

        if col == 'Unnamed: 0':
            continue

        value = df_row[col].values[0]

        if not isinstance(value,str) and np.isnan(value):
            value = None

        data[col] = value

    return data