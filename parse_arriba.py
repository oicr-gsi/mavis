#!/usr/bin/env python3

import csv
import argparse

parser = argparse.ArgumentParser(description='output file')
parser.add_argument("path", type=str, help='input file (Arriba results)')
parser.add_argument('-o', '--output', help='output file path', required=True)
args = parser.parse_args()


def parse_arriba(file):
    fusion_list = []
    with open(file, newline='') as fusions:
        fusion_reader = csv.DictReader(fusions, delimiter='\t')
        for row in fusion_reader:
            fusion_list.append(dict(row))
    fusions.close()
    return fusion_list


# ============ Here we start creating an array of output lines ===================================
lines = ["\t".join(['#break1_chromosome', 'break1_position_start', 'break1_position_end',
                    'break2_chromosome', 'break2_position_start', 'break2_position_end', 'tools']) + "\n"]
# ============ Parse lines from Arriba calls =====================================================
arriba_lines = parse_arriba(args.path)

for line in arriba_lines:
    br1 = line.get('breakpoint1').split(":")
    br2 = line.get('breakpoint2').split(":")
    lines.append("\t".join([br1[0], br1[1], br1[1], br2[0], br2[1], br2[1], "arriba"]) + "\n")

with open(args.output, mode='+w') as out:
    out.writelines(lines)
out.close()
