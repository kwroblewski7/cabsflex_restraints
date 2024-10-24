import argparse
import os
import json
import re
import numpy as np
from biopandas.pdb import PandasPdb

def make_restraints(pdb_file, plddt_file, ss_file, working_directory, output_file):
    check_work_dir(working_directory)
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"File '{pdb_file}' not found")
    pdb_df = parse_pdb(pdb_file)
    if plddt_file == 'bf':
        pdb_df = pdb_df.rename(columns={'b_factor': 'plddt'})
    else:
        if not os.path.isfile(plddt_file):
            raise FileNotFoundError(f"File '{plddt_file}' not found")
        pdb_df = pdb_df.drop(columns=['b_factor'])
        plddts = parse_plddt(plddt_file)
        pdb_df['plddt'] = plddts

    if not os.path.isfile(ss_file):
        raise FileNotFoundError(f"File '{ss_file}' not found")
    sec = parse_ss(ss_file)
    pdb_df['ss'] = pdb_df.apply(lambda row: sec[(str(row['residue_number']), row['chain_id'])], axis=1)

    pdb_df['category'] = pdb_df.apply(category_methods, axis=1)

    generate_restraints(pdb_df, output_file)

def check_work_dir(working_directory):
    working_directory = working_directory.rstrip('/')
    parent_directory = os.path.dirname(working_directory)
    if not parent_directory:
        parent_directory = '.'
    if os.path.isdir(parent_directory):
        os.makedirs(working_directory, exist_ok=True)
    else:
        raise FileNotFoundError(f"Directory '{parent_directory}' not found")
    os.chdir(working_directory)

def parse_pdb(pdb_file):
    pdb = PandasPdb().read_pdb(pdb_file)
    df = pdb.df['ATOM'][pdb.df['ATOM']['atom_name'] == 'CA']
    return df[['residue_number', 'chain_id' ,'x_coord', 'y_coord', 'z_coord', 'b_factor']].reset_index(drop=True)

def parse_plddt(plddt_file):
    with open(plddt_file, 'r') as file:
        plddt_values = json.load(file)
    plddts = plddt_values['plddt']
    return plddts

def parse_ss(ss_file):
    with open(ss_file, 'r') as file:
        sec = {}
        sec_line = '^([0-9 ]{5}) ([0-9 ]{4}.)([A-Z ]) ([A-Z])  ([HBEGITSP ])(.*)$'

        for line in file.read().split('\n'):
            matched_line = re.match(sec_line, line)
            if matched_line:
                residue_number = str(matched_line.group(2).strip())
                chain = str(matched_line.group(3))
                if matched_line.group(5) in 'HGIP':
                    ss = 'H'
                elif matched_line.group(5) in 'BE':
                    ss = 'E'
                elif matched_line.group(5) in 'T':
                    ss = 'T'
                else:
                    ss = 'C'
                sec[(residue_number, chain)] = ss
    return sec

def category_methods(row):

    category = 0
    if row['plddt'] <= 50:
        pass
    elif row['plddt'] <= 70:
        category += 1
    elif row['plddt'] <= 90:
        category += 2
    else:
        category += 3

    if row['ss'] == 'C':
        category -= 1
    elif row['ss'] == 'T':
        pass
    else:
        category += 1

    if category < 0:
        category = 0
    elif category > 3:
        category = 3

    return category

def generate_restraints(
        pdb_df, output_file, distance = (3.8, 11.5), distance_width = 1.0,
        sequence_gap = 3, restraint_strength = (3.5, 0.5)
):
    with open(output_file, 'w') as output:
        for i in range(len(pdb_df)):
            for j in range(i + sequence_gap, len(pdb_df)):
                res1 = pdb_df.iloc[i]
                res2 = pdb_df.iloc[j]
                if res1['chain_id'] == res2['chain_id']:
                    restraint_distance = _calculate_pairwise_distance(res1[['x_coord', 'y_coord', 'z_coord']],
                                                                      res2[['x_coord', 'y_coord', 'z_coord']])
                    if distance[0] <= restraint_distance <= distance[1]:
                        category_sum = res1['category'] + res2['category']
                        if category_sum < 4:
                            weight = 0.0
                        elif category_sum == 4:
                            weight = 0.5
                        else:
                            weight = 1.0
                        strength_min = restraint_strength[0]  * weight
                        strength_max = restraint_strength[1]  * weight
                        if weight:
                            restraint_string = (
                                "{residue1}:{chain_id} {residue2}:{chain_id} {distance:.2f} {distance_width:.2f} "
                                "{min_weights:.2f} {max_weights:.2f}\n"
                            ).format(
                                residue1=res1['residue_number'],
                                residue2=res2['residue_number'],
                                chain_id=res1['chain_id'],
                                distance=restraint_distance,
                                distance_width=distance_width,
                                min_weights=strength_min,
                                max_weights=strength_max
                            )
                            output.write(restraint_string)

def _calculate_pairwise_distance(res1_coords, res2_coords):
    return np.linalg.norm(np.array(res1_coords) - np.array(res2_coords))

def parser():
    parser_instance = argparse.ArgumentParser(description='Generate restraints data from pLDDT scores and secondary structure')
    parser_instance.add_argument('-i', '--pdb_file', type=str, required=True, help='Path to PDB file')
    parser_instance.add_argument('-p', '--plddt_file', type=str, required=True, help="Can be either:\n - Path to pLDDT score file\n - 'bf' to use B-factor values as pLDDT scores")
    parser_instance.add_argument('-s','--ss_file', type=str, required=True, help='Path to secondary structure file')
    parser_instance.add_argument('-w', '--work_dir', type=str, default='.', help='Working directory. All files paths should be relative to this directory')
    parser_instance.add_argument('-o', '--output_file', type=str, default='ca_restraints_file.txt', help='Output file name')
    arguments = parser_instance.parse_args()
    return arguments


if __name__ == '__main__':
    args = parser()
    make_restraints(args.pdb_file, args.plddt_file, args.ss_file, args.work_dir, args.output_file)