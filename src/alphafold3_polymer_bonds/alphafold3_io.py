
import argparse, collections, copy, gzip, importlib, importlib.metadata, importlib.resources, json, os, os.path, re, string, sys
from copy import deepcopy
from pathlib import Path
from pprint import pprint
from typing import Dict, List, Tuple, Any, TypeAlias

import humanfriendly

JSON: TypeAlias = dict[str, "JSON"] | list["JSON"] | str | int | float | bool | None #https://github.com/python/typing/issues/182#issuecomment-1320974824

def _open_r(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')

def _encode_indices_arrays(js):
    #https://github.com/google-deepmind/alphafold3/blob/v3.0.1/src/alphafold3/common/folding_input.py#L1294-L1302
    return re.sub(
        r'("(?:queryIndices|templateIndices)": \[)([\s\n\d,]+)(\],?)',
        lambda mtch: mtch[1] + re.sub(r'\n\s+', ' ', mtch[2].strip()) + mtch[3],
        js,
    )

def read_input_json(path):
    """Read json while preserving order of keys from file"""
    with _open_r(path) as fh:
        return json.load(fh, object_pairs_hook=collections.OrderedDict)

def print_json(js, max_size=500):
    """Print json without long MSA strings"""
    js_ = copy.deepcopy(js)
    for sequence in js_['sequences']:
        seq_type = next(iter(sequence.keys()))
        if seq_type in {'protein', 'dna', 'rna'}:
            for k, v in sequence[seq_type].items():
                if len(v) > max_size:
                    sequence[seq_type][k] = f'<{humanfriendly.format_size(len(v))} string>'
    print(json.dumps(js_, indent=2))

def write_input_json(js, path):
    """Write json aiming to match AF3; if path contains {}, replaces with name from js"""
    js_str = _encode_indices_arrays(json.dumps(js, indent=2))
    if '{}' in path:
        path = path.format(js['name'])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(js_str)

def santise_name(s):
    """AF3 job names are lower case, numeric, -._
    https://github.com/google-deepmind/alphafold3/blob/v3.0.1/src/alphafold3/common/folding_input.py#L857-L861
    """
    def is_allowed(c):
        return c.islower() or c.isnumeric() or c in set('-._')
    return ''.join(filter(is_allowed, s.strip().lower().replace(' ', '_')))

def count_tokens(path):
    """Count tokens
    TODO: proteins only, no nucleic acids, no ligands, no PTMs...
    """
    sequences = read_input_json(path)['sequences']
    n_tokens = 0
    for seq in sequences:
        if 'protein' in seq:
            n_chains = len(seq['protein']['id'])
            seq_len = len(seq['protein']['sequence'])
            n_tokens += n_chains * seq_len
    return n_tokens

def multimer_json(*monomers):
    js = copy.deepcopy(monomers[0])
    js['name'] = '_'.join([monomer['name'] for monomer in monomers])
    for monomer in monomers[1:]:
        js['sequences'].append(copy.deepcopy(monomer['sequences'][0]))
    for monomer, chain_id in zip(js['sequences'], string.ascii_uppercase):
        monomer['protein']['id'] = [chain_id]
    return js

def read_summary_confidences(path, name):
    js = read_input_json(os.path.join(path, name, f'{name}_summary_confidences.json'))
    return js
