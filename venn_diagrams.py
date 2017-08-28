#! /usr/bin/env python3

import pandas as pd
from matplotlib_venn import venn3_unweighted

def prepare_kraken(kraken):

    kraken['SampleName'] = kraken['Sample'].str.split('_').str[0]+'_'+kraken['Sample'].str.split('_').str[1]

    kraken['Type'] = kraken['SampleName'].str.split('_').str[1]

    kraken['Type'] = kraken['Type'].str.replace('2','')

    kraken[kraken['Kraken.ID'] == 374840]

    # Remove phiX and microviridae reads

    kraken = kraken.drop(kraken['Kraken.ID'] == 374840)
    kraken = kraken.drop(kraken['Kraken.ID'] == 10842)
    kraken = kraken.drop(kraken['Kraken.ID'] == 10841)

# Remove

def categoryDiffKraken(category):
    full = kraken.loc[(kraken['Level.of.Classification'] == category) & (kraken['Type'] =="D1"), "Name"]
    half = kraken.loc[(kraken['Level.of.Classification'] == category) & (kraken['Type'] == "D0.5"), "Name"]
    quar = kraken.loc[(kraken['Level.of.Classification'] == category) & (kraken['Type'] == "D0
      .25"), "Name"]
    setfull = set(full)
    sethalf = set(half)
    setquar = set(quar)
    return(setfull, sethalf, setquar)

def categoryDiffAMR(category):
    full = amrResults[category][amrResults['Sample_type'] == "D1"]
    half = amrResults[category][amrResults['Sample_type'] == "D0.5"]
    quar = amrResults[category][amrResults['Sample_type'] == "D0.25"]
    seqtk = amrResults[category][amrResults['Sample_type'] == "D0.25_seqtk"]
    setfull = set(full)
    sethalf = set(half)
    setquar = set(quar)
    setseqtk = set(seqtk)
    return(setfull, sethalf, setquar, setseqtk)

allCatsSets = {}

for c in allCats:
    allCatsSets[c] = categoryDiffKraken(c)

familyFullSet = allCatsSets['F'][0]
familyHalfSet = allCatsSets['F'][1]
familyQuarSet = allCatsSets['F'][2]
familySets = [familyFullSet, familyHalfSet, familyQuarSet]

classFullSet = allCatsSets['C'][0]
classHalfSet = allCatsSets['C'][1]
classQuarSet = allCatsSets['C'][2]
classSets = [classFullSet, classHalfSet, classQuarSet]

orderFullSet = allCatsSets['O'][0]
orderHalfSet = allCatsSets['O'][1]
orderQuarSet = allCatsSets['O'][2]
orderSets = [orderFullSet, orderHalfSet, orderQuarSet]

phylumFullSet = allCatsSets['P'][0]
phylumHalfSet = allCatsSets['P'][1]
phylumQuarSet = allCatsSets['P'][2]
phylumsets = [phylumfullset, phylumhalfset, phylumquarset]


genusFullSet = allCatsSets['g'][0]
genusHalfSet = allCatsSets['g'][1]
genusQuarSet = allCatsSets['g'][2]
genusSets = [genusFullSet, genusHalfSet, genusQuarSet]


v3Phylum = venn3_unweighted(phylumSets, ('D1', 'D0.5', 'D0.25'))
plt.title('Phylum')
plt.savefig('phylumVennUpdated2.png')
plt.clf()
plt.cla()

v3Order = venn3_unweighted(orderSets, ('D1', 'D0.5', 'D0.25'))
plt.title('Order')
plt.savefig('orderVennUpdated2.png')
plt.clf()
plt.cla()

v3Family = venn3_unweighted(familySets, ('D1', 'D0.5', 'D0.25'))
plt.title('Family')
plt.savefig('familyVennUpdated2.png')
plt.clf()
plt.cla()

v3Class = venn3_unweighted(classSets, ('D1', 'D0.5', 'D0.25'))
plt.title('Class')
plt.savefig('amrclassVennUpdated2.png')
plt.clf()
plt.cla()

v3Genus = venn3_unweighted(genusSets, ('D1', 'D0.5', 'D0.25'))
plt.title('Genus')
plt.savefig('genusVennUpdated2.png')
plt.clf()
plt.cla()

# Building AMR Venn diagrams

allAMRCats = ['Class', 'Mechanism', 'Group', 'Name']

allAMRCatsSets = {}

for c in allAMRCats:
    allAMRCatsSets[c] = categoryDiffAMR(c)

classFullSet = allAMRCatsSets['Class'][0]
classHalfSet = allAMRCatsSets['Class'][1]
classQuarSet = allAMRCatsSets['Class'][2]
classSeqtkSet = allAMRCatsSets['Class'][3]
classSets = [classFullSet, classHalfSet, classQuarSet, classSeqtkSet]

mechFullSet = allAMRCatsSets['Mechanism'][0]
mechHalfSet = allAMRCatsSets['Mechanism'][1]
mechQuarSet = allAMRCatsSets['Mechanism'][2]
mechSeqtkSet = allAMRCatsSets['Mechanism'][3]
mechSets = [mechFullSet, mechHalfSet, mechQuarSet, mechSeqtkSet]

groupFullSet = allAMRCatsSets['Group'][0]
groupHalfSet = allAMRCatsSets['Group'][1]
groupQuarSet = allAMRCatsSets['Group'][2]
groupSeqtkSet = allAMRCatsSets['Group'][3]
groupSets = [groupFullSet, groupHalfSet, groupQuarSet, groupSeqtkSet]

geneFullSet = allAMRCatsSets['Name'][0]
geneHalfSet = allAMRCatsSets['Name'][1]
geneQuarSet = allAMRCatsSets['Name'][2]
geneSeqtkSet = allAMRCatsSets['Name'][3]
geneSets = [geneFullSet, geneHalfSet, geneQuarSet, geneSeqtkSet]

v3class = venn3_unweighted(classSets[0:2], ('D1', 'D0.5', 'D0.25'))
plt.title('Class')
plt.savefig('amrClassVenn.png')
plt.clf()
plt.cla()

v3mech = venn3_unweighted(mechSets[0:2], ('D1', 'D0.5', 'D0.25'))
plt.title('Mech')
plt.savefig('amrMechVenn.png')
plt.clf()
plt.cla()

v3group = venn3_unweighted(groupSets[0:2], ('D1', 'D0.5', 'D0.25'))
plt.title('Group')
plt.savefig('amrGroupVenn.png')
plt.clf()
plt.cla()

v3gene = venn3_unweighted(geneSets[0:2], ('D1', 'D0.5', 'D0.25'))
plt.title('Gene')
plt.savefig('amrGeneVenn.png')
plt.clf()
plt.cla()


def setOperations():
    """Generate a file with all the results of the set operations
    Ideally return a table with the results"""

geneAllIntersects = geneFullSet.intersection(geneHalfSet).intersection(geneQuarSet)
geneFullDiffHalf = geneFullSet.difference(geneHalfSet)
geneFullDiffQuar = geneFullSet.difference(geneQuarSet)
geneFullInterHalfDiffQuar = geneFullSet.intersection(geneHalfSet).difference(geneQuarSet)
geneHalfDiffFull = geneHalfSet.difference(geneFullSet)
geneHalfDiffQuar = geneHalfSet.difference(geneQuarSet)
geneHalfInterQuarDiffFull = geneHalfSet.intersection(geneFullSet).difference(geneFullSet)
geneQuarDiffFull = geneQuarSet.difference(geneFullSet)
geneQuarDiffHalf = geneQuarSet.difference(geneHalfSet)
geneQuarInterFullDiffHalf = geneQuarSet.intersection(geneFullSet).difference(geneHalfSet)


classAllIntersects = classFullSet.intersection(classHalfSet).intersection(classQuarSet)
classFullDiffHalf = classFullSet.difference(classHalfSet)
classFullDiffQuar = classFullSet.difference(classQuarSet)
classFullInterHalfDiffQuar = classFullSet.intersection(classHalfSet).difference(classQuarSet)
classHalfDiffFull = classHalfSet.difference(classFullSet)
classHalfDiffQuar = classHalfSet.difference(classQuarSet)
classHalfInterQuarDiffFull = classHalfSet.intersection(classFullSet).difference(classFullSet)
classQuarDiffFull = classQuarSet.difference(classFullSet)
classQuarDiffHalf = classQuarSet.difference(classHalfSet)
classQuarInterFullDiffHalf = classQuarSet.intersection(classFullSet).difference(classHalfSet)


groupAllIntersects = groupFullSet.intersection(groupHalfSet).intersection(groupQuarSet)
groupFullDiffHalf = groupFullSet.difference(groupHalfSet)
groupFullDiffQuar = groupFullSet.difference(groupQuarSet)
groupFullInterHalfDiffQuar = groupFullSet.intersection(groupHalfSet).difference(groupQuarSet)
groupHalfDiffFull = groupHalfSet.difference(groupFullSet)
groupHalfDiffQuar = groupHalfSet.difference(groupQuarSet)
groupHalfInterQuarDiffFull = groupHalfSet.intersection(groupFullSet).difference(groupFullSet)
groupQuarDiffFull = groupQuarSet.difference(groupFullSet)
groupQuarDiffHalf = groupQuarSet.difference(groupHalfSet)
groupQuarInterFullDiffHalf = groupQuarSet.intersection(groupFullSet).difference(groupHalfSet)


mechAllIntersects = mechFullSet.intersection(mechHalfSet).intersection(mechQuarSet)
mechFullDiffHalf = mechFullSet.difference(mechHalfSet)
mechFullDiffQuar = mechFullSet.difference(mechQuarSet)
mechFullInterHalfDiffQuar = mechFullSet.intersection(mechHalfSet).difference(mechQuarSet)
mechHalfDiffFull = mechHalfSet.difference(mechFullSet)
mechHalfDiffQuar = mechHalfSet.difference(mechQuarSet)
mechHalfInterQuarDiffFull = mechHalfSet.intersection(mechFullSet).difference(mechFullSet)
mechQuarDiffFull = mechQuarSet.difference(mechFullSet)
mechQuarDiffHalf = mechQuarSet.difference(mechHalfSet)
mechQuarInterFullDiffHalf = mechQuarSet.intersection(mechFullSet).difference(mechHalfSet)

geneSets = [geneAllIntersects,
geneFullDiffHalf,
geneFullDiffQuar,
geneFullInterHalfDiffQuar,
geneHalfDiffFull,
geneHalfDiffQuar,
geneHalfInterQuarDiffFull,
geneQuarDiffFull,
geneQuarDiffHalf,
geneQuarInterFullDiffHalf]

classSets = [classAllIntersects,
classFullDiffHalf,
classFullDiffQuar,
classFullInterHalfDiffQuar,
classHalfDiffFull,
classHalfDiffQuar,
classHalfInterQuarDiffFull,
classQuarDiffFull,
classQuarDiffHalf,
classQuarInterFullDiffHalf]

mechSets = [mechAllIntersects,
mechFullDiffHalf,
mechFullDiffQuar,
mechFullInterHalfDiffQuar,
mechHalfDiffFull,
mechHalfDiffQuar,
mechHalfInterQuarDiffFull,
mechQuarDiffFull,
mechQuarDiffHalf,
mechQuarInterFullDiffHalf]

groupSets = [groupAllIntersects,
groupFullDiffHalf,
groupFullDiffQuar,
groupFullInterHalfDiffQuar,
groupHalfDiffFull,
groupHalfDiffQuar,
groupHalfInterQuarDiffFull,
groupQuarDiffFull,
groupQuarDiffHalf,
groupQuarInterFullDiffHalf]


geneSetLengths = [len(g) for g in geneSets]
classSetLengths = [len(c) for c in classSets]
mechSetLengths = [len(m) for m in mechSets]
groupSetLengths = [len(g) for g in groupSets]

geneLists = [list(g) for g in geneSets]
classLists = [list(c) for c in classSets]
mechLists = [list(m) for m in mechSets]
groupLists = [list(g) for g in groupSets]

amrGeneSetOperationsTuple = zip(setOperationKeys, geneSetLengths, geneLists)
amrClassSetOperationsTuple = zip(setOperationKeys, classSetLengths, classLists)
amrMechSetOperationsTuple = zip(setOperationKeys, mechSetLengths, mechLists)
amrGroupSetOperationsTuple = zip(setOperationKeys, groupSetLengths, groupLists)

amrGeneUniqueDF = pd.DataFrame.from_records(amrGeneSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
amrClassUniqueDF = pd.DataFrame.from_records(amrClassSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
amrMechUniqueDF = pd.DataFrame.from_records(amrMechSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
amrGroupUniqueDF = pd.DataFrame.from_records(amrGroupSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])

amrGeneUniqueDF.to_csv('amrGeneSetOperations.csv', index = False)
amrClassUniqueDF.to_csv('amrClassSetOperations.csv', index = False)
amrMechUniqueDF.to_csv('amrMechSetOperations.csv', index = False)
amrGroupUniqueDF.to_csv('amrGroupSetOperations.csv', index = False)


