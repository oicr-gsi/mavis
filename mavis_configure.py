#!/usr/bin/python

from optparse import OptionParser

"""MAVIS configurator
                
This script prepares a bash script that may be used to configure
MAVIS, annotation tool for structural variants (SV). It will accept
multiple inputs and form the body of the script with three types 
of instructions:

 --library   this is to specify analysis type using designated inputs
 --convert   for SV callers
 --assign
        
This tool accepts comma separated values .
            
"""
USAGE = "mavis_configure.py -bam [comma-separated bams] -lib [comma-separated lib designs] -svdata [comma-separated files with SV variants] -workflows [comma-separated names of workflows]"

libtypes = {'WT': "transcriptome", 'MR': "transcriptome", 'WG': "genome"}
assign_arrays = {'WT': [], 'MR': [], 'WG': []}
wfMappings = {'StructuralVariation': 'delly', 'Delly': 'delly', 'StarFusion': 'starfusion', 'Manta': 'manta'}
svMappings = {'delly': ["WG"], 'starfusion': ["WT","MR"], 'manta': ['WG']}

parser = OptionParser()
parser.add_option("-b", "--bam", dest="bam",
                  help="input bams")
parser.add_option("-l", "--lib", dest="lib",
                  help="input library names")
parser.add_option("-s", "--svdata", dest="svdata",
                  help="input variants")
parser.add_option("-w", "--workflow", dest="wf",
                  help="workflow names")
parser.add_option("-d", "--donor", dest="donor",
                  help="donor id")
parser.add_option("-c", "--config", dest="conf",
                  help="config file")

# Small subroutine to check the validity of inputs

def validate_options(parser: OptionParser) -> OptionParser:
    (opts, arguments) = parser.parse_args()
    """Make sure we have bams, sv calls and matching number of metadata pieces"""
    if opts.bam is None:
        print("We need input bam files")
        exit(1)
    if opts.lib is None:
        print("We need library types")
        exit(1)
    if opts.svdata is None:
        print("We need files with SV data")
        exit(1)
    if opts.wf is None:
        print("We need workflow names")
        exit(1)
    if opts.conf is None:
        print("We need a file to write to")
        exit(1)
    if opts.donor is None:
        print("We need a donor id")
        exit(1)
    bams = opts.bam.split(",")
    libs = opts.lib.split(",")
    svdata = opts.svdata.split(",")
    wfs = opts.wf.split(",")
    if len(bams) != len(libs) or len(svdata) != len(wfs):
        print("Lengths of data array don not match")
        exit(1)
    return parser


# 1. Make sure we have all the right options

ok_parser = validate_options(parser)

# 2. Parse inputs and print lines into stdout
(options, args) = parser.parse_args()

bams = options.bam.split(",")
libs = options.lib.split(",")
svdata = options.svdata.split(",")
wfs = options.wf.split(",")

# Line Arrays:
library_lines = []
convert_lines = []
assign_lines = []

for b in range(len(bams)):
    flag = ('False' if libs[b] == 'WG' else 'True')
    library_lines.append("--library " + libs[b] + "." + options.donor + " " + libtypes[libs[b]] + " diseased " + flag + " " + bams[b])

for s in range(len(svdata)):
    for w in wfMappings.keys():
        if w in wfs[s]:
            convert_lines.append("--convert " + wfMappings[w] + " " + svdata[s] + " " + wfMappings[w])
            for library_type in svMappings[wfMappings[w]]:
                assign_arrays[library_type].append(wfMappings[w])

# Parse again bam files, make assign lines
for b in range(len(bams)):
    if len(assign_arrays[libs[b]]) > 0:
        separator = " "
        tools = separator.join(assign_arrays[libs[b]])
        assign_lines.append("--assign " + libs[b] + "." + options.donor + " " + tools)

# 3. print to stdout
print("#!/bin/bash" + "\n")
print('mavis config ' + ' \\')
for lib_line in library_lines:
    print(lib_line + ' \\')
for convert_line in convert_lines:
    print(convert_line + ' \\')
for assign_line in assign_lines:
    print(assign_line + ' \\')

print("--write " + options.conf)

