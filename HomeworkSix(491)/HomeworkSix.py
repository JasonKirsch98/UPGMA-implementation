
### Parse fasta files, store with hashtable

def read_fasta(input_file):
    taxon, seq = None, []
    for line in input_file:
        line = line.rstrip()
        if line.startswith(">"):
            if taxon: yield (taxon, ''.join(seq))
            taxon, seq = line[1:], []
        else:
            seq.append(line)
    if taxon: yield (taxon, ''.join(seq))




### compute Hamming distance matrix
#compare each sequence to every other, +1 if character match, 0 if not.
def computeDistance(seq1, seq2):
    dist = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            dist+=1
    return dist


def getUpgmaTableIndex(tax1, tax2):
    if tax1 > tax2:
        return tax1 + '|' + tax2
    return tax2+'|'+tax1


def main():
    suspect_sequence = ''
    with open('query_file') as input_file:
        for name, seq in read_fasta(input_file):
            suspect_sequence = seq

    # Read database file
    database = {}
    with open('database_file') as input_file:
        for name, seq in read_fasta(input_file):
            database[name] = seq
    database['SUSPECT'] = suspect_sequence

    # Create UPGMA data table
    upgmaTable = {}
    activeSequences = {}
    numActiveSequences = 0
    for keyA, valueA in database.iteritems():
        activeSequences[keyA] = 1
        numActiveSequences += 1
        for keyB, valueB in database.iteritems():
            if keyA > keyB:
                distance = computeDistance(valueA, valueB)
                index = getUpgmaTableIndex(keyA, keyB)
                upgmaTable[index] = distance
    while numActiveSequences > 2:
        # Find smallest distance in upgma table
        shortestDistance = 10000000
        shortestDistanceKey = ''
        for key,value in upgmaTable.iteritems():
            if value < shortestDistance:
                shortestDistance = value
                shortestDistanceKey = key
        commaIndex = shortestDistanceKey.index('|')

        # Get both keys
        keyA = shortestDistanceKey[0:commaIndex]
        keyB = shortestDistanceKey[commaIndex+1:]
        activeSequences[keyA] = 0
        activeSequences[keyB] = 0
        combinedKey = '(' +keyA+','+keyB+')'
        upgmaTable.pop(shortestDistanceKey, None)

        # Remove original keys from table, add combined keys to table
        newActiveSequences = []
        for activeKey, value in activeSequences.iteritems():
            if value == 1:
                upgmaIndexA = getUpgmaTableIndex(activeKey, keyA)
                upgmaIndexB = getUpgmaTableIndex(activeKey, keyB)
                upgmaValA = upgmaTable.pop(upgmaIndexA, None)
                upgmaValB = upgmaTable.pop(upgmaIndexB, None)

                averageDistance = (upgmaValA + upgmaValB) / 2
                combinedUpgmaIndex = getUpgmaTableIndex(combinedKey, activeKey)
                upgmaTable[combinedUpgmaIndex] = averageDistance
        activeSequences[combinedKey] = 1
        numActiveSequences -= 1

    # Get the resulting tree expression

    # One remaining key in the UPGMA table
    lastKey = ''
    for key, value in upgmaTable.iteritems():
        lastKey = key

    # Split remaining two nodes
    commaIndex = shortestDistanceKey.index('|')
    keyA = shortestDistanceKey[0:commaIndex]
    keyB = shortestDistanceKey[commaIndex + 1:]

    resultingTreeExpression = '(' + keyA + ',' + keyB + ");"

    print(resultingTreeExpression)


if __name__ == '__main__':
    main()