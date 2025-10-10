
import argparse, collections, collections.abc, copy, hashlib, gzip, json, os, os.path, re, string, subprocess, sys
from copy import deepcopy
from pathlib import Path
from pprint import pprint
from typing import Dict, List, Tuple, Any, TypeAlias

import Bio.PDB, Bio.PDB.mmcifio, Bio.PDB.Polypeptide

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

def _sequence_hash(seq):
    return hashlib.sha1(seq.encode()).hexdigest()

def init_input_json(*seqs):
    def _get_seq(id, seq):
        return collections.OrderedDict([('protein', collections.OrderedDict([('id', id),('sequence', seq)]))])
    js = collections.OrderedDict([
        ('dialect', 'alphafold3'),
        ('version', 2),
        ('name', 'name'),
        ('sequences', []),
        ('modelSeeds', [1]),
        ('bondedAtomPairs', None),
        ('userCCD', None)])
    for seq, chain_id in zip(seqs, string.ascii_uppercase):
        js['sequences'].append(_get_seq(chain_id, seq))
    return js

def read_input_json(path):
    """Read json while preserving order of keys from file"""
    with _open_r(path) as fh:
        return json.load(fh, object_pairs_hook=collections.OrderedDict)

def print_input_json(js, max_size=500):
    """Print (part of) json without long MSA strings"""
    def iter_(js):
        if isinstance(js, str) or isinstance(js, int) or isinstance(js, list):
            return js
        for k, v in js.items():
            if k in {'templates', 'unpairedMsa', 'pairedMsa'} and len(v) > max_size:
                js[k] = f'<{humanfriendly.format_size(len(v))} string>'
            elif isinstance(v, collections.abc.Mapping):
                js[k] = iter_(v)
            elif isinstance(v, list):
                for i in range(len(v)):
                    v[i] = iter_(v[i])
        return js
    print(json.dumps(iter_(js), indent=2))

def write_input_json(js, path):
    """Write json aiming to match AF3; if path contains {}, replaces with name from js"""
    # Infer path from name attribute
    if '{}' in path:
        path = path.format(js['name'])

    # Infer name attribute from path
    basename = os.path.basename(path).rstrip('.json')
    if not(basename.startswith(js['name'])):
        js['name'] = basename
        print(f'Inferring name attribute {js["name"]} from path - {path}')

    if os.path.dirname(path) != '':
        os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        js_str = _encode_indices_arrays(json.dumps(js, indent=2))
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

def read_summary_confidences(path, name):
    js = read_input_json(os.path.join(path, name, f'{name}_summary_confidences.json'))
    return js

def get_colabfold_msa(seq, dir='/tmp/_get_colabfold_msa'):
    name = _sequence_hash(seq)
    path_input = f'{dir}/input/{name}.fasta'
    #print(path_input)
    os.makedirs(os.path.dirname(path_input), exist_ok=True)
    with open(path_input, 'w') as f:
        f.write(f'>{name}\n{seq}')

    path_output = f'{dir}/output/{name}.json'
    if not os.path.isfile(path_output):
        path_output_dir = f'{dir}/output'
        os.makedirs(path_output_dir, exist_ok=True)
        cmd = f'MPLBACKEND=AGG; source /colabfold_venv/bin/activate; colabfold_batch --msa-only --af3-json {path_input} {path_output_dir}'
        print(cmd)
        r = subprocess.run(cmd, capture_output=True, shell=True, executable='/bin/bash')
        assert r.returncode == 0

    return read_input_json(path_output)

def colab_data_pipeline(js):
    for seq in js['sequences']:
        if 'protein' in seq.keys():
            seq_msa = get_colabfold_msa(seq['protein']['sequence'])['sequences'][0]['protein']
            for field in ['templates', 'unpairedMsa', 'pairedMsa']:
                seq['protein'][field] = seq_msa[field]
    return js

def get_structure(path, only_first=True):
    """Attempt to read a structure using Bio.PDB while transparently handling compression/PDB-vs-CIF"""
    # Check file format - pdb1 is used by rcsb bioassemblies...
    if path.endswith('.pdb') or path.endswith('.pdb1') or path.endswith('.pdb.gz') or path.endswith('.pdb1.gz'):
        parser = Bio.PDB.PDBParser(QUIET=True)
    else:
        parser = Bio.PDB.MMCIFParser(QUIET=True)

    # Check for .gz compresssion
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as fh:
            struct = parser.get_structure(path, fh)
    else:
        struct = parser.get_structure(path, path)

    return struct[0] if only_first else struct

def input_json_sequence_init(type, id, sequence):
    return collections.OrderedDict([
        (type, collections.OrderedDict([
            ('id', id),
            ('sequence', sequence)
        ]))
    ])

_dna_dict = {
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DT': 'T',
    'TGP': 'G',
}

def input_json_sequence_from_chain(chain):
    # Filter out waters & print het-entities
    resnames = []
    for residue in chain.get_residues():
        if residue.get_id()[0] != 'W':
            resnames.append(residue.resname)
            if residue.get_id()[0] != ' ':
                print(chain.id, residue.resname, residue.get_id())

    # Check if entities are a subset with <=
    if set(resnames) <= {'DA', 'DC', 'DG', 'DT', 'TGP', 'TTG', 'TTC'}:
        type = 'dna'
        sequence = ''.join([ _dna_dict.get(resname, '') for resname in resnames ])
    elif set(resnames) <= {'A', 'C', 'G', 'U'}:
        type = 'rna'
        sequence = ''.join(resnames)
    else:
        type = 'protein'
        sequence = ''.join([ Bio.PDB.Polypeptide.protein_letters_3to1.get(resname, '') for resname in resnames ])

    return input_json_sequence_init(
        type=type,
        id=chain.id,
        sequence=sequence,
    )

def input_json_from_rcsb(path):
    struct = get_structure(path)
    js = init_input_json()
    js['name'] = os.path.basename(path).removesuffix('.gz').removesuffix('.cif').removesuffix('.pdb')
    for chain in Bio.PDB.Selection.unfold_entities(entity_list=struct, target_level='C'):
        js['sequences'].append(input_json_sequence_from_chain(chain))
    return js
