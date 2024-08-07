from PDNPR.PDNPR import pdnpr

Mutation_AA = ['G132A', 'R130P', 'Y68N', 'C124S', 'C136S', 'R130L', 'R130Q', 'Y155C']
md_file = [
    'md-out/high-G132A/fit.xtc',
    'md-out/high-R130P/fit.xtc',
    'md-out/high-Y68N/fit.xtc',
    'md-out/or-C124S/fit.xtc',
    'md-out/or-C136S/fit.xtc',
    'md-out/or-R130L/fit.xtc',
    'md-out/or-R130Q/fit.xtc',
    'md-out/or-Y155C/fit.xtc'
]
pdb_file = [
    'md-out/high-G132A/P.pdb',
    'md-out/high-R130P/P.pdb',
    'md-out/high-Y68N/P.pdb',
    'md-out/or-C124S/P.pdb',
    'md-out/or-C136S/P.pdb',
    'md-out/or-R130L/P.pdb',
    'md-out/or-R130Q/P.pdb',
    'md-out/or-Y155C/P.pdb'
]

step = 10
end_AA = 127
edge_cutoff = 0.5


for i in range(len(Mutation_AA)):
    print("===============================================")
    print(f"============ Handling Mutation {Mutation_AA[i]} ==========")
    print("===============================================")
    pdnpr(step, Mutation_AA[i], end_AA, edge_cutoff, md_file[i], pdb_file[i])
