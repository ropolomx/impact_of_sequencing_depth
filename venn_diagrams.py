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
geneAllInter = geneFullSet.intersection(geneHalfSet).intersection(geneQuarSet)


