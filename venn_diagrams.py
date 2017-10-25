#! /usr/bin/env python3

import pandas as pd
from pandas import ExcelWriter
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted

def read_kraken(report):

    kraken = pd.read_csv(report)

    return kraken

def prepare_kraken(kraken):

    kraken["SampleName"] = kraken['Sample'].str.split('_').str[0] + '_'

    kraken['Sample'].str.split('_').str[1]

    kraken['Type'] = kraken['SampleName'].str.split('_').str[1]

    kraken['Type'] = kraken['Type'].str.replace('2','')

    # Remove phiX and microviridae reads
    # Update this section according to the needs and discoveries of kraken datasets

    removal = [374840, # phiX
               10842, # Microvirinae
               10841 # Microviridae
               ]

    kraken = kraken[~kraken['Kraken.ID'].isin(removal)]

    return kraken


# Remove

def categoryDiffKraken(category):
    full = krakenConcat.loc[(krakenConcat['TaxRank'] == category) & (krakenConcat['Sample_Type'] =="F"), "Name"]
    half = krakenConcat.loc[(krakenConcat['TaxRank'] == category) & (krakenConcat['Sample_Type'] == "H"), "Name"]
    quar = krakenConcat.loc[(krakenConcat['TaxRank'] == category) & (krakenConcat['Sample_Type'] == "QD"), "Name"]
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

allTaxa = ['P', 'C', 'F', 'O', 'G', 'S']

allTaxaSets = {}

for t in allTaxa:
    allTaxaSets[t] = categoryDiffKraken(t)

familyFullSet = allTaxaSets['F'][0]
familyHalfSet = allTaxaSets['F'][1]
familyQuarSet = allTaxaSets['F'][2]
familySets = [familyFullSet, familyHalfSet, familyQuarSet]

classFullSet = allTaxaSets['C'][0]
classHalfSet = allTaxaSets['C'][1]
classQuarSet = allTaxaSets['C'][2]
classSets = [classFullSet, classHalfSet, classQuarSet]

orderFullSet = allTaxaSets['O'][0]
orderHalfSet = allTaxaSets['O'][1]
orderQuarSet = allTaxaSets['O'][2]
orderSets = [orderFullSet, orderHalfSet, orderQuarSet]

phylumFullSet = allTaxaSets['P'][0]
phylumHalfSet = allTaxaSets['P'][1]
phylumQuarSet = allTaxaSets['P'][2]
phylumSets = [phylumFullSet, phylumHalfSet, phylumQuarSet]


genusFullSet = allTaxaSets['G'][0]
genusHalfSet = allTaxaSets['G'][1]
genusQuarSet = allTaxaSets['G'][2]
genusSets = [genusFullSet, genusHalfSet, genusQuarSet]


speciesFullSet = allTaxaSets['S'][0]
speciesHalfSet = allTaxaSets['S'][1]
speciesQuarSet = allTaxaSets['S'][2]
speciesSets = [speciesFullSet, speciesHalfSet, speciesQuarSet]

# Add option to plot friendlier Venn diagrams for color-blind audiences

v3Phylum = venn3_unweighted(phylumSets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Phylum.set_labels:
    text.set_fontsize(16)
for text in v3Phylum.subset_labels:
    text.set_fontsize(18)
plt.title('Phylum', fontsize=20)
plt.savefig('krakenPhylumVennCB.png')
plt.clf()
plt.cla()

v3Order = venn3_unweighted(orderSets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Order.set_labels:
    text.set_fontsize(16)
for text in v3Order.subset_labels:
    text.set_fontsize(18)
plt.title('Order', fontsize=20)
plt.savefig('krakenOrderVennCB.png')
plt.clf()
plt.cla()

v3Family = venn3_unweighted(familySets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Family.set_labels:
    text.set_fontsize(16)
for text in v3Family.subset_labels:
    text.set_fontsize(18)
plt.title('Family', fontsize=20)
plt.savefig('krakenFamilyVennCB.png')
plt.clf()
plt.cla()

v3Class = venn3_unweighted(classSets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Class.set_labels:
    text.set_fontsize(16)
for text in v3Class.subset_labels:
    text.set_fontsize(18)
plt.title('Class', fontsize=20)
plt.savefig('krakenClassVennCB.png')
plt.clf()
plt.cla()

v3Genus = venn3_unweighted(genusSets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Genus.set_labels:
    text.set_fontsize(16)
for text in v3Genus.subset_labels:
    text.set_fontsize(18)
plt.title('Genus',fontsize=20)
plt.savefig('krakenGenusVennCB.png')
plt.clf()
plt.cla()

v3Species = venn3_unweighted(speciesSets, ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3Species.set_labels:
    text.set_fontsize(16)
for text in v3Species.subset_labels:
    text.set_fontsize(18)
plt.title('Species', fontsize=20)
plt.savefig('krakenSpeciesVennCB.png')
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

v3class = venn3_unweighted(classSets[0:3], ('D1', 'D0.5', 'D0.25'))
for text in v3class.set_labels:
    text.set_fontsize(16)
for text in v3class.subset_labels:
    text.set_fontsize(18)
plt.title('Class', fontsize=20)
plt.savefig('amrClassVenn.png')
plt.clf()
plt.cla()

v3class = venn3_unweighted(classSets[0:3], ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3class.set_labels:
    text.set_fontsize(16)
for text in v3class.subset_labels:
    text.set_fontsize(18)
plt.title('Class', fontsize=20)
plt.savefig('amrClassVennCB.png')
plt.clf()
plt.cla()

v3mech = venn3_unweighted(mechSets[0:3], ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3mech.set_labels:
    text.set_fontsize(16)
for text in v3mech.subset_labels:
    text.set_fontsize(18)
plt.title('Mechanism', fontsize=20)
plt.savefig('amrMechVennCB.png')
plt.clf()
plt.cla()

v3group = venn3_unweighted(groupSets[0:3], ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3group.set_labels:
    text.set_fontsize(16)
for text in v3group.subset_labels:
    text.set_fontsize(18)
plt.title('Group', fontsize=20)
plt.savefig('amrGroupVennCB.png')
plt.clf()
plt.cla()

v3gene = venn3_unweighted(geneSets[0:3], ('D1', 'D0.5', 'D0.25'), set_colors=('y', 'b', 'r'))
for text in v3gene.set_labels:
    text.set_fontsize(16)
for text in v3gene.subset_labels:
    text.set_fontsize(18)
plt.title('Gene', fontsize=20)
plt.savefig('amrGeneVennCB.png')
plt.clf()
plt.cla()


setOperationKeys = ['All_intersections',
        'Full_vs_Half',
        'Full_vs_Quar',
        'Full_and_Half_vs_Quar',
        'Half_vs_Full',
        'Half_vs_Quar',
        'Half_and_Quar_vs_Full',
        'Quar_vs_Full',
        'Quar_vs_Half',
        'Quar_and_Full_vs_Half']


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

allClassComps = []
allMechanismComps = []
allGroupComps = []
allGeneComps = []

for c in classLists:
    allClassComps.append(amrResults[amrResults['Class'].isin(c)])

for m in mechLists:
    allMechanismComps.append(amrResults[amrResults['Mechanism'].isin(m)])

for g in groupLists:
    allGroupComps.append(amrResults[amrResults['Group'].isin(g)])

for g in geneLists:
    allGeneComps.append(amrResults[amrResults['Name'].isin(g)])

classSlices = zip(setOperationKeys, allClassComps)

mechSlices = zip(setOperationKeys, allMechanismComps)

groupSlices = zip(setOperationKeys, allGroupComps)

geneSlices = zip(setOperationKeys, allGeneComps)


classWriter = ExcelWriter('classComparisons.xlsx')
mechWriter = ExcelWriter('mechComparisons.xlsx')
groupWriter = ExcelWriter('groupComparisons.xlsx')
geneWriter = ExcelWriter('geneComparisons.xlsx')

for n, df in classSlices:
    df.to_excel(classWriter, sheet_name=n, index=False)
classWriter.save()

for n, df in mechSlices:
    df.to_excel(mechWriter, sheet_name=n, index=False)
mechWriter.save()

for n, df in groupSlices:
    df.to_excel(groupWriter, sheet_name=n, index=False)
groupWriter.save()

for n, df in geneSlices:
    df.to_excel(geneWriter, sheet_name=n, index=False)
geneWriter.save()

# Kraken data analysis


def krakenSetOperations():
    """Generate a file with all the results of the set operations
    Ideally return a table with the results"""

speciesAllIntersects = speciesFullSet.intersection(speciesHalfSet).intersection(speciesQuarSet)
speciesFullDiffHalf = speciesFullSet.difference(speciesHalfSet)
speciesFullDiffQuar = speciesFullSet.difference(speciesQuarSet)
speciesFullInterHalfDiffQuar = speciesFullSet.intersection(speciesHalfSet).difference(speciesQuarSet)
speciesHalfDiffFull = speciesHalfSet.difference(speciesFullSet)
speciesHalfDiffQuar = speciesHalfSet.difference(speciesQuarSet)
speciesHalfInterQuarDiffFull = speciesHalfSet.intersection(speciesFullSet).difference(speciesFullSet)
speciesQuarDiffFull = speciesQuarSet.difference(speciesFullSet)
speciesQuarDiffHalf = speciesQuarSet.difference(speciesHalfSet)
speciesQuarInterFullDiffHalf = speciesQuarSet.intersection(speciesFullSet).difference(speciesHalfSet)


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


genusAllIntersects = genusFullSet.intersection(genusHalfSet).intersection(genusQuarSet)
genusFullDiffHalf = genusFullSet.difference(genusHalfSet)
genusFullDiffQuar = genusFullSet.difference(genusQuarSet)
genusFullInterHalfDiffQuar = genusFullSet.intersection(genusHalfSet).difference(genusQuarSet)
genusHalfDiffFull = genusHalfSet.difference(genusFullSet)
genusHalfDiffQuar = genusHalfSet.difference(genusQuarSet)
genusHalfInterQuarDiffFull = genusHalfSet.intersection(genusFullSet).difference(genusFullSet)
genusQuarDiffFull = genusQuarSet.difference(genusFullSet)
genusQuarDiffHalf = genusQuarSet.difference(genusHalfSet)
genusQuarInterFullDiffHalf = genusQuarSet.intersection(genusFullSet).difference(genusHalfSet)


familyAllIntersects = familyFullSet.intersection(familyHalfSet).intersection(familyQuarSet)
familyFullDiffHalf = familyFullSet.difference(familyHalfSet)
familyFullDiffQuar = familyFullSet.difference(familyQuarSet)
familyFullInterHalfDiffQuar = familyFullSet.intersection(familyHalfSet).difference(familyQuarSet)
familyHalfDiffFull = familyHalfSet.difference(familyFullSet)
familyHalfDiffQuar = familyHalfSet.difference(familyQuarSet)
familyHalfInterQuarDiffFull = familyHalfSet.intersection(familyFullSet).difference(familyFullSet)
familyQuarDiffFull = familyQuarSet.difference(familyFullSet)
familyQuarDiffHalf = familyQuarSet.difference(familyHalfSet)
familyQuarInterFullDiffHalf = familyQuarSet.intersection(familyFullSet).difference(familyHalfSet)

orderAllIntersects = orderFullSet.intersection(orderHalfSet).intersection(orderQuarSet)
orderFullDiffHalf = orderFullSet.difference(orderHalfSet)
orderFullDiffQuar = orderFullSet.difference(orderQuarSet)
orderFullInterHalfDiffQuar = orderFullSet.intersection(orderHalfSet).difference(orderQuarSet)
orderHalfDiffFull = orderHalfSet.difference(orderFullSet)
orderHalfDiffQuar = orderHalfSet.difference(orderQuarSet)
orderHalfInterQuarDiffFull = orderHalfSet.intersection(orderFullSet).difference(orderFullSet)
orderQuarDiffFull = orderQuarSet.difference(orderFullSet)
orderQuarDiffHalf = orderQuarSet.difference(orderHalfSet)
orderQuarInterFullDiffHalf = orderQuarSet.intersection(orderFullSet).difference(orderHalfSet)

phylumAllIntersects = phylumFullSet.intersection(phylumHalfSet).intersection(phylumQuarSet)
phylumFullDiffHalf = phylumFullSet.difference(phylumHalfSet)
phylumFullDiffQuar = phylumFullSet.difference(phylumQuarSet)
phylumFullInterHalfDiffQuar = phylumFullSet.intersection(phylumHalfSet).difference(phylumQuarSet)
phylumHalfDiffFull = phylumHalfSet.difference(phylumFullSet)
phylumHalfDiffQuar = phylumHalfSet.difference(phylumQuarSet)
phylumHalfInterQuarDiffFull = phylumHalfSet.intersection(phylumFullSet).difference(phylumFullSet)
phylumQuarDiffFull = phylumQuarSet.difference(phylumFullSet)
phylumQuarDiffHalf = phylumQuarSet.difference(phylumHalfSet)
phylumQuarInterFullDiffHalf = phylumQuarSet.intersection(phylumFullSet).difference(phylumHalfSet)

phylumSets = [phylumAllIntersects,
phylumFullDiffHalf,
phylumFullDiffQuar,
phylumFullInterHalfDiffQuar,
phylumHalfDiffFull,
phylumHalfDiffQuar,
phylumHalfInterQuarDiffFull,
phylumQuarDiffFull,
phylumQuarDiffHalf,
phylumQuarInterFullDiffHalf]

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

orderSets = [orderAllIntersects,
orderFullDiffHalf,
orderFullDiffQuar,
orderFullInterHalfDiffQuar,
orderHalfDiffFull,
orderHalfDiffQuar,
orderHalfInterQuarDiffFull,
orderQuarDiffFull,
orderQuarDiffHalf,
orderQuarInterFullDiffHalf]

familySets = [familyAllIntersects,
familyFullDiffHalf,
familyFullDiffQuar,
familyFullInterHalfDiffQuar,
familyHalfDiffFull,
familyHalfDiffQuar,
familyHalfInterQuarDiffFull,
familyQuarDiffFull,
familyQuarDiffHalf,
familyQuarInterFullDiffHalf]

genusSets = [genusAllIntersects,
genusFullDiffHalf,
genusFullDiffQuar,
genusFullInterHalfDiffQuar,
genusHalfDiffFull,
genusHalfDiffQuar,
genusHalfInterQuarDiffFull,
genusQuarDiffFull,
genusQuarDiffHalf,
genusQuarInterFullDiffHalf]


speciesSets = [speciesAllIntersects,
speciesFullDiffHalf,
speciesFullDiffQuar,
speciesFullInterHalfDiffQuar,
speciesHalfDiffFull,
speciesHalfDiffQuar,
speciesHalfInterQuarDiffFull,
speciesQuarDiffFull,
speciesQuarDiffHalf,
speciesQuarInterFullDiffHalf]

phylumSetLengths = [len(p) for p in phylumSets]
classSetLengths = [len(c) for c in classSets]
orderSetLengths = [len(o) for o in orderSets]
familySetLengths = [len(f) for f in familySets]
genusSetLengths = [len(g) for g in genusSets]
speciesSetLengths = [len(s) for s in speciesSets]

phylumSetLists = [list(p) for p in phylumSets]
classSetLists = [list(c) for c in classSets]
orderSetLists = [list(o) for o in orderSets]
familySetLists = [list(f) for f in familySets]
genusSetLists = [list(g) for g in genusSets]
speciesSetLists = [list(s) for s in speciesSets]

setOperationKeys = ['All_intersections',
        'Full_vs_Half',
        'Full_vs_Quar',
        'Full_and_Half_vs_Quar',
        'Half_vs_Full',
        'Half_vs_Quar',
        'Half_and_Quar_vs_Full',
        'Quar_vs_Full',
        'Quar_vs_Half',
        'Quar_and_Full_vs_Half']


phylumSetOperationsTuple = zip(setOperationKeys, phylumSetLengths, phylumSetLists)
classSetOperationsTuple = zip(setOperationKeys, classSetLengths, classSetLists)
orderSetOperationsTuple = zip(setOperationKeys, orderSetLengths, orderSetLists)
familySetOperationsTuple = zip(setOperationKeys, familySetLengths, familySetLists)
genusSetOperationsTuple = zip(setOperationKeys, genusSetLengths, genusSetLists)
speciesSetOperationsTuple = zip(setOperationKeys, speciesSetLengths, speciesSetLists)

phylumUniqueDF = pd.DataFrame.from_records(phylumSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
classUniqueDF = pd.DataFrame.from_records(classSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
orderUniqueDF = pd.DataFrame.from_records(orderSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
familyUniqueDF = pd.DataFrame.from_records(familySetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
genusUniqueDF = pd.DataFrame.from_records(genusSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])
speciesUniqueDF = pd.DataFrame.from_records(speciesSetOperationsTuple, columns=['Comparison', 'Size', 'Members'])

phylumUniqueDF.to_csv('phylumSetOperations.csv', index = False)
classUniqueDF.to_csv('classSetOperations.csv', index = False)
orderUniqueDF.to_csv('orderSetOperations.csv', index = False)
familyUniqueDF.to_csv('familySetOperations.csv', index = False)
genusUniqueDF.to_csv('genusSetOperations.csv', index = False)
speciesUniqueDF.to_csv('speciesSetOperations.csv', index = False)

allPhylumComps = []
allClassComps = []
allOrderComps = []
allFamilyComps = []
allGenusComps = []
allSpeciesComps = []

for p in phylumSetLists:
    allPhylumComps.append(krakenConcat[krakenConcat['Name'].isin(p)])

for c in classSetLists:
    allClassComps.append(krakenConcat[krakenConcat['Name'].isin(c)])

for o in orderSetLists:
    allOrderComps.append(krakenConcat[krakenConcat['Name'].isin(o)])

for f in familySetLists:
    allFamilyComps.append(krakenConcat[krakenConcat['Name'].isin(f)])

for g in genusSetLists:
    allGenusComps.append(krakenConcat[krakenConcat['Name'].isin(g)])

for s in speciesSetLists:
    allSpeciesComps.append(krakenConcat[krakenConcat['Name'].isin(s)])

phylumSlices = zip(setOperationKeys, allPhylumComps)

classSlices = zip(setOperationKeys, allClassComps)

orderSlices = zip(setOperationKeys, allOrderComps)

familySlices = zip(setOperationKeys, allFamilyComps)

genusSlices = zip(setOperationKeys, allGenusComps)

speciesSlices = zip(setOperationKeys, allSpeciesComps)

phylumWriter = ExcelWriter('phylumComparisons.xlsx')
classWriter = ExcelWriter('classComparisons.xlsx')
orderWriter = ExcelWriter('orderComparisons.xlsx')
familyWriter = ExcelWriter('familyComparisons.xlsx')
genusWriter = ExcelWriter('genusComparisons.xlsx')
speciesWriter = ExcelWriter('speciesComparisons.xlsx')

for n, df in phylumSlices:
    df.to_excel(phylumWriter, sheet_name=n, index=False)
phylumWriter.save()

for n, df in classSlices:
    df.to_excel(classWriter, sheet_name=n, index=False)
classWriter.save()

for n, df in orderSlices:
    df.to_excel(orderWriter, sheet_name=n, index=False)
orderWriter.save()

for n, df in familySlices:
    df.to_excel(familyWriter, sheet_name=n, index=False)
familyWriter.save()

for n, df in genusSlices:
    df.to_excel(genusWriter, sheet_name=n, index=False)
genusWriter.save()

for n, df in speciesSlices:
    df.to_excel(speciesWriter, sheet_name=n, index=False)
speciesWriter.save()
