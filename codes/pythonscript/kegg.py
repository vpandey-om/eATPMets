# Standard library packages
import io
import os

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

# Import Pandas, so we can use dataframes
import pandas as pd


def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

result = REST.kegg_list("pathway", "eco").read()
df=to_df(result)

result = REST.kegg_get("pbe00740").read()
print(result)

result = REST.kegg_link("compound", "map00740").read()
compund_df=to_df(result)

import pdb; pdb.set_trace()
