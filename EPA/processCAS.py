import os, logging
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
PandasTools.RenderImagesInAllDataFrames()

from tqdm import tqdm_notebook
tqdm_notebook().pandas()

from EPA import comptox_lookup, comptox_link

# Parse arguments
import argparse
parser = argparse.ArgumentParser(description='Resolve CAS numbers into smiles and calculate molecular descriptors.')
parser.add_argument('-f', '--filename', type=argparse.FileType('r'), help='File containing CAS numbers')
parser.add_argument('-c', '--column', type=int, help='Column in the input file containing CAS numbers')
args = parser.parse_args()

# Set up a logger object...
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('[%(asctime)s %(name)s %(levelname)s] %(message)s', datefmt='%d/%m/%y %H:%M:%S'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# Load CAS numbers
#filename = "C:\\Users\\bet\\Documents\\Inditex\\Compound inventory\\Annotated compounds\\REACH\\AnnexVI_inventory\\20170516\\annex_vi_clp_table_en.split.TallAndSkinny.tsv"
with open(args.filename, encoding='utf8') as infile:
    cas_numbers = set([x.split('\t')[args.column-1] for x in infile])

# Lookup CAS numbers
records = [y for y in (comptox_lookup(x) for x in tqdm_notebook(cas_numbers)) if y is not None]

# Load information into a dataframe
comptox_df = pd.DataFrame(records)
# Check for missing names or structures...
for x in ('name', 'inchi', 'inchikey', 'smiles'): assert (comptox_df[x].isnull() | (comptox_df[x] == '')).sum() == 0
comptox_df['report'] = comptox_df['cas'].apply(comptox_link)   # Add link to CompTox dashboard...

# Add RDKit molecule and descriptors...
comptox_df['mol'] = [Chem.MolFromMolBlock(x) for x in comptox_df['conntab']]
comptox_df['MW'] = [Descriptors.MolWt(x) if not x == None else None for x in comptox_df['mol']]
comptox_df['logP'] = [Descriptors.MolLogP(x) if not x == None else None for x in comptox_df['mol']]
comptox_df['TPSA'] = [Descriptors.TPSA(x) if not x == None else None for x in comptox_df['mol']]

# Dump dataframe into output file
#outFile = "C:\\Users\\bet\\Documents\\Inditex\\Compound inventory\\Annotated compounds\\REACH\\AnnexVI_inventory\\20170516\\annex_vi_clp_table_en.split.TallAndSkinny.smi_desc.tsv"
bn, ext = os.path.splitext(args.filename)
outFile = bn + '.smiles_desc.tsv'
outF = open(outFile, 'w')
comptox_df.drop(['inchi', 'conntab', 'synonyms', 'mol'], axis=1, inplace=True)
comptox_df.to_csv(outFile, sep='\t', float_format='%.2f', index=False)
outF.close()