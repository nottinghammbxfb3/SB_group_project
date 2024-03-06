def mutation_counter(ref_protein, ref_name, alt_protein, alt_name, out_file):
    '''
    Get amount, position and mutations between reference and alternative protein strings
    In: 2 EQUAL LENGTH protein sequences, names of sequence origin and output file (.txt)
    Out: .txt file containing mutation data
    '''
    count=0
    mutations=[]

    # Identify mutated amino acids and append information to list

    for i in range(len(ref_protein)):
        if ref_protein[i] != alt_protein[i]:
            count += 1
            mutations.append([i + 1, ref_protein[i], alt_protein[i]])

    with open(out_file, 'w') as f:
        f.write(f'Number of mutations: {count}\n\n')
        for pos, ref, alt in mutations:
            f.write(f'Pos: {pos}\n')
            f.write(f'{ref_name}: {ref}\n')
            f.write(f'{alt_name}: {alt}\n\n')

# Example input:
'''
mutation_counter(
    'MSVVVFSPMYPSSHWCLDELDKLFAINKRMNHTMLPIFYKVDPSHVRKPNVHFKDDFDRHVSDGRFNEEKMERWRRAMIDF_',
    'Ref',
    'MSVVVFSPTYPSSHWCLDELDKLFAINKRMNHTMLPIFYKVNPSHVRKPNVHFKDDFDRHVSDGRFNEEKMERWRRAMIDF_',
    'Alt',
    'test_seq_align_out.txt')
'''
