from Bio.PDB.DSSP import make_dssp_dict
import re

name = "test"
sequence = "KECV"
fastas = [[name, sequence]]

print("############ ASA ####################")

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    header.append('ASA.F' + str(p))
encodings.append(header)

with open('Seq_99.spXout') as f:
    records = f.readlines()[1:]

for i in fastas:

    code = [name]
    proteinSeq = ''
    asaValue = []
    for line in records:
        array = line.strip().split() if line.strip() != '' else None # a single line of spxout file
        proteinSeq = proteinSeq + array[1] # the amino acids from spxout file, aggregated to a peptide sequence
        asaValue.append(array[10]) # the actual asa value
    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        for p in range(pos, pos + len(sequence)):
            code.append(asaValue[p])
    encodings.append(code)

print(encodings)

###

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    header.append('ASA.F' + str(p))
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
for k, v in dssp[0].items():
    # print(f"{k}: {v[2]}")
    code.append(v[2])

encodings.append(code)

print(encodings)

print("############ TA ####################")

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    header.append('TA.F' + str(p) + '.phi')
    header.append('TA.F' + str(p) + '.psi')
encodings.append(header)

for i in fastas:
    name, sequence = i[0], i[1]
    code = [name]

    with open('Seq_99.spXout') as f:
        records = f.readlines()[1:]

    proteinSeq = ''
    asaValue = []
    for line in records:
        array = line.strip().split() if line.strip() != '' else None
        proteinSeq = proteinSeq + array[1]
        asaValue.append(array[3:5])
    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        for p in range(pos, pos+len(sequence)):
            code.append(asaValue[p][0])
            code.append(asaValue[p][1])
    encodings.append(code)

print(encodings)

###

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    header.append('TA.F' + str(p) + '.phi')
    header.append('TA.F' + str(p) + '.psi')
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
for k, v in dssp[0].items():
    code.append(v[3])
    code.append(v[4])

encodings.append(code)

print(encodings)

print("############ SSEC ####################")


def calculateSSE(pos, end, SSE):
    newValues = SSE[pos:end]
    return [newValues.count('H')/(end-pos), newValues.count('E')/(end-pos), newValues.count('C')/(end-pos)]


encodings = []
header = ['#', 'H', 'E', 'C']
encodings.append(header)

for i in fastas:
    name, sequence = i[0], i[1]
    code = [name]

    with open('Seq_99.spXout') as f:
        records = f.readlines()[1:]

    proteinSeq = ''
    SSE = []
    for line in records:
        array = line.strip().split() if line.rstrip() != '' else None
        proteinSeq = proteinSeq + array[1]
        SSE.append(array[2])

    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        code = code + calculateSSE(pos, pos+len(sequence), SSE)
    encodings.append(code)

print(encodings)

###

mappings = {
    "H": "H", "G": "H", "I": "H",
    "S": "C", "T": "C", "C": "C", "-": "C",
    "E": "E", "B": "E"
}

encodings = []
header = ['#', 'H', 'E', 'C']
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
H, E, C = 0, 0, 0
for k, v in dssp[0].items():
    if mappings[v[1]] == "H":
        H += 1
    elif mappings[v[1]] == "E":
        E += 1
    elif mappings[v[1]] == "C":
        C += 1

seq_len = len(dssp[0].keys())
code += [H/seq_len, E/seq_len, C/seq_len]

encodings.append(code)

print(encodings)

print("############ SSEB ####################")

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    for ss in ('H', 'E', 'C'):
            header.append('Pos'+str(p)+'.'+ss)
encodings.append(header)

for i in fastas:
    name, sequence = i[0], i[1]
    code = [name]

    with open('Seq_99.spXout') as f:
        records = f.readlines()[1:]

    proteinSeq = ''
    SSE = []
    myDict = {'H':[0, 0, 1], 'E':[0, 1, 0], 'C':[1, 0, 0]}
    for line in records:
        array = line.strip().split() if line.rstrip() != '' else None
        proteinSeq = proteinSeq + array[1]
        SSE.append(array[2])

    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        for p in range(pos, pos + len(sequence)):
            code = code + myDict[SSE[p]]
    encodings.append(code)

print(encodings)

###

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    for ss in ('H', 'E', 'C'):
            header.append('Pos'+str(p)+'.'+ss)
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
myDict = {'H':[1, 0, 0], 'E':[0, 1, 0], 'C':[0, 0, 1]}
for k, v in dssp[0].items():
    if mappings[v[1]] == "H":
        code += myDict["H"]
    elif mappings[v[1]] == "E":
        code += myDict["E"]
    elif mappings[v[1]] == "C":
        code += myDict["C"]

encodings.append(code)

print(encodings)

print("############ DISORDER B ####################")

encodings = []
header = ['#']
for p in range(1, 2*len(fastas[0][1])+1):
    header.append('disorderB.F' + str(p))

encodings.append(header)
for i in fastas:
    name, sequence = i[0], i[1]
    code = [name]

    with open("Seq_99.dis") as f:
        records = f.readlines()

    tag = 0
    for i in range(len(records)):
        if re.search('^-------', records[i]):
            tag = i
            break
    records = records[tag+1:-1]

    proteinSeq = ''
    disValue = []
    myDict = {'D':[0, 1], 'O':[1, 0]}
    for line in records:
        array = line.rstrip().split() if line.rstrip() != '' else None
        key = array[3] if array[3] == 'D' else 'O'
        proteinSeq = proteinSeq + array[1]
        disValue.append(key)

    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        for p in range(pos, pos+len(sequence)):
            code = code + myDict[disValue[p]]
    encodings.append(code)

print(encodings)

###

encodings = []
header = ['#']
for p in range(1, 2*len(fastas[0][1])+1):
    header.append('disorderB.F' + str(p))

encodings = []
header = ['#']
for p in range(1, len(fastas[0][1])+1):
    for ss in ('H', 'E', 'C'):
            header.append('Pos'+str(p)+'.'+ss)
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
myDict = {'D': [1, 0], 'O': [0, 1]}
for k, v in dssp[0].items():
    if v[1] == "-":
        code += myDict["D"]
    else:
        code += myDict["O"]

encodings.append(code)

print(encodings)

print("############ DISORDER C ####################")


def calculateDicorderContent(pos, endPos, disValue):
    newValues = disValue[pos: endPos]
    return [newValues.count('D')/(endPos - pos), newValues.count('O')/(endPos-pos)]


encodings = []
header = ['#', 'disorder-content', 'order-content']
encodings.append(header)

for i in fastas:
    name, sequence = i[0], i[1]
    code = [name]

    with open("Seq_99.dis") as f:
        records = f.readlines()

    tag = 0
    for i in range(len(records)):
        if re.search('^-------', records[i]):
            tag = i
            break
    records = records[tag + 1:-1]

    proteinSeq = ''
    disValue = []
    for line in records:
        array = line.rstrip().split() if line.rstrip() != '' else None
        key = array[3] if array[3] == 'D' else 'O'
        proteinSeq = proteinSeq + array[1]
        disValue.append(key)

    pos = proteinSeq.find(sequence)
    if pos == -1:
        print('Warning: could not find the peptide in proteins.\n\n')
    else:
        code = code + calculateDicorderContent(pos, pos + len(sequence), disValue)
    encodings.append(code)

print(encodings)

###

encodings = []
header = ['#', 'disorder-content', 'order-content']
encodings.append(header)

code = [name]
dssp = make_dssp_dict("Seq_4.dssp")
D, O = 0, 0
for k, v in dssp[0].items():
    if v[1] == "-":
        D += 1
    else:
        O += 1

seq_len = len(dssp[0].keys())
code += [D/seq_len, O/seq_len]

encodings.append(code)

print(encodings)
