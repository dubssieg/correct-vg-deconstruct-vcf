from sys import argv
from os.path import exists
from os import remove,system

gfa_file_path,vcf_file_path,output_vcf,reference_name,alternative_name = argv[1:]

def revcomp(string: str, compl: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}) -> str:
    try:
        return ''.join([compl[s] for s in string.upper()][::-1])
    except IndexError as exc:
        raise IndexError(
            "Complementarity does not include all chars in sequence.") from exc
    
if not exists("rs-pancat-paths"):
    system("wget https://github.com/dubssieg/rs-pancat-paths/releases/download/0.1.3/rs-pancat-paths")
    system("chmod +x rs-pancat-paths")

index_gfa:str = f"{gfa_file_path}.index"
masked_gfa:str = f"{gfa_file_path}.mask.gfa"
masked_vcf:str = f"{vcf_file_path}.mask.vcf"

system(f"./rs-pancat-paths {gfa_file_path} index > {index_gfa}")

paths_to_mask:list[str] = list()
with open(index_gfa,'r') as reader:
    for line in reader:
        if not line.startswith('#'):
            if not line.split('\t')[0] in [reference_name,alternative_name]:
                paths_to_mask.append(f"-M {line.split('\t')[0]}")

system(f"./rs-pancat-paths {gfa_file_path} mask {' '.join(paths_to_mask)} > {masked_gfa}")
system(f"vg deconstruct -a {masked_gfa} -p {reference_name} > {masked_vcf}")
    
reference_gfa_set:set = set()
alternative_gfa_set:set = set()
sequences_dict:dict = dict()

with open(gfa_file_path,'r',encoding='utf-8') as gfa_reader:
    for line in gfa_reader:
        if line.startswith('S'):
            sequences_dict[line.split('\t')[1]] = line.strip().split('\t')[2]
        if line.startswith('P') and reference_name in line.split('\t')[1]:
            for x in line.split('\t')[2].split(','):
                reference_gfa_set.add(int(x[:-1]))
        if line.startswith('P') and alternative_name in line.split('\t')[1]:
            for x in line.split('\t')[2].split(','):
                alternative_gfa_set.add(int(x[:-1]))

private_alternative_set:set = alternative_gfa_set - reference_gfa_set

alternative_vcf_set:set = set()
headers:list = list()
alt_pos:int = 0
incorrects_formattings:list[str] = list()

with open(vcf_file_path,'r',encoding='utf-8') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            headers = line.split('\t')
            alt_pos = headers.index(alternative_name)
        else:
            # on check si le genome alt est décrit en génotype par un entier 1...n
            try:
                has_variant:int = int(line.split('\t')[alt_pos])
            except:
                incorrects_formattings.append(line.split('\t')[alt_pos])
                has_variant:int = -1
            # si c'est 0 ou -1, c'est que ce n'est pas correctement décrit comme variant
            if has_variant > 0:
                # On doit aller chercher la bonne position du variant dans le champ AT (pas besoin de conserver le sens)
                alternative_vcf_set.update([int(x) for x in line.split('\t')[7].split(';')[3].split(',')[has_variant].replace('<','>')[1:].split('>')])

resulting_set:set = private_alternative_set - alternative_vcf_set

alternative_gfa_list:list = list()
alternative_gfa_strs:list = list()

with open(gfa_file_path,'r',encoding='utf-8') as gfa_reader:
    for line in gfa_reader:
        if line.startswith('P') and alternative_name in line.split('\t')[1]:
            for x in line.split('\t')[2].split(','):
                alternative_gfa_list.append(int(x[:-1]))
                alternative_gfa_strs.append(x)

intersect:set = alternative_gfa_set.intersection(reference_gfa_set)

count_rescuables_nodes:int = 0
node_pairs:set[tuple] = set()
chains_list:list = list()

for node_id in resulting_set:
    idx:int = alternative_gfa_list.index(node_id)
    bubble_branch:list = [alternative_gfa_strs[idx]]
    # We assume we have a single occurence, we search the first index
    cutoff:int = min(int(len(alternative_gfa_list)/10),idx,len(alternative_gfa_list)-idx) -1
    x,y = 0,0
    dist = 1
    while not x and cutoff >= dist:
        bubble_branch = [alternative_gfa_strs[idx-dist]] + bubble_branch
        if alternative_gfa_list[idx-dist] in intersect:
            x = alternative_gfa_list[idx-dist]
        dist+=1
    dist = 1
    while not y and cutoff >= dist:
        bubble_branch = bubble_branch + [alternative_gfa_strs[idx+dist]]
        if alternative_gfa_list[idx+dist] in intersect:
            y = alternative_gfa_list[idx+dist]
        dist+=1
    if x and y:
        count_rescuables_nodes+=1
        node_pairs.add((x,y))
        chains_list.append(bubble_branch)
print(f"Graph contains at least {count_rescuables_nodes} rescuable nodes ({round(count_rescuables_nodes/len(resulting_set)*100,2)}%).")

reference_gfa_list:list = list()
reference_gfa_strs:list = list()

with open(gfa_file_path,'r',encoding='utf-8') as gfa_reader:
    for line in gfa_reader:
        if line.startswith('P') and reference_name in line.split('\t')[1]:
            for x in line.split('\t')[2].split(','):
                reference_gfa_list.append(int(x[:-1]))
                reference_gfa_strs.append(x)

bubbles_dict = dict()
bubbles_positions = dict()
bubbles_positions_alt = dict()

for a,b in node_pairs:
    idx_a,idx_b = reference_gfa_list.index(a),reference_gfa_list.index(b)
    bubbles_dict[(a,b)] = reference_gfa_strs[min(idx_a,idx_b):max(idx_a+1,idx_b+1)]
    bubbles_positions[(a,b)] = sum([len(sequences_dict[x[:-1]]) for x in reference_gfa_strs[0:min(idx_a,idx_b)]])
    bubbles_positions_alt[(a,b)] = sum([len(sequences_dict[x[:-1]]) for x in alternative_gfa_strs[0:min(idx_a,idx_b)]])

reference_chains:list[list[int]] = list()

with open(vcf_file_path,'r',encoding='utf-8') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('##') or line.startswith('#'):
            continue
        else:
            # on récupère la chaîne de référence
            reference_chains.append([int(x) for x in line.split('\t')[7].split(';')[3].split(',')[0].replace('<','>')[4:].split('>')])

masked_bubbles:list = list()

with open(masked_vcf,'r') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('#'):
            continue
        else:
            masked_bubbles.append(line.split('\t')[2])

if exists(output_vcf):
    remove(output_vcf)

pairlist = list()

counter:int = 0
for bubble in chains_list:
    bubble_name = ('>' if bubble[0][-1] == '+' else '<')+bubble[0][:-1]+('>' if bubble[-1][-1] == '+' else '<')+bubble[-1][:-1]
    bubble_anchors:set[int] = {int(bubble[0][:-1]),int(bubble[1][:-1])}
    intricate = [len(bubble_anchors.intersection(x[1:-1])) for x in reference_chains]
    if sum(intricate) == 0:
        status = '0'
    else:
        intricate = sorted([x for x in intricate if x != 0])
        if intricate[0] == 1 and intricate[-1] == 2:
            status = '1|2'
        else:
            status = str(intricate[0])

    if bubble_name in pairlist:
        continue
    else:
        counter+=1
        pairlist.append(bubble_name)
        print(
            reference_name+
            '\t'+
            str(bubbles_positions[(int(bubble[0][:-1]),int(bubble[-1][:-1]))])+
            '\t'+
            str(bubbles_positions_alt[(int(bubble[0][:-1]),int(bubble[-1][:-1]))])+
            '\t'+
            bubble_name+
            '\t'+
            ''.join([sequences_dict[idd[:-1]] if idd[-1] == '+' else revcomp(sequences_dict[idd[:-1]]) for idd in bubbles_dict[(int(bubble[0][:-1]),int(bubble[-1][:-1]))]][1:-1])+
            '\t'+
            ''.join([sequences_dict[idd[:-1]] if idd[-1] == '+' else revcomp(sequences_dict[idd[:-1]]) for idd in bubble][1:-1])+
            '\t'+'60'+
            '\t'+'.'+
            "\t"+
            "AC="+";"+
            "AF="+";"+
            "AN="+";"+
            "AT="+
            ''.join([(">" if idd[-1] == '+' else '<')+str(idd)[:-1] for idd in bubbles_dict[(int(bubble[0][:-1]),int(bubble[-1][:-1]))]])+','+''.join([(">" if idd[-1] == '+' else '<')+str(idd)[:-1] for idd in bubble])+";"+
            "NS="+";"+
            "LV="+";"+
            "RS="+str(int(bubble_name in masked_bubbles))+";"+
            "IT="+status
            ,file=open(output_vcf,'a+')
        )

print(f"Sucessfully rescued {counter} bubbles!")

# Clean temp files
if exists(masked_vcf):
    remove(masked_vcf)

if exists(masked_gfa):
    remove(masked_gfa)

if exists(index_gfa):
    remove(index_gfa)
