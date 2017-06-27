#! /usr/bin/env python3

import pandas as pd

kraken['SampleName'] = kraken['Sample'].str.split('_').str[0]+'_'+kraken['Sample'].str.split('_').str[1]

kraken['Type'] = kraken['SampleName'].str.split('_').str[1]

kraken['Type'] = kraken['Type'].str.replace('2','')

kraken[kraken['Kraken.ID'] == 374840]

kraken = kraken.drop(kraken['Kraken.ID'] == 10842)
kraken = kraken.drop(kraken['Kraken.ID'] == 10841)

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
    setfull = set(full)
    sethalf = set(half)
    setquar = set(quar)
    return(setfull, sethalf, setquar)

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



geneFullSet = allAMRCatsSets['Name'][0]
geneHalfSet = allAMRCatsSets['Name'][1]
geneQuarSet = allAMRCatsSets['Name'][2]
geneSets = [geneFullSet, geneHalfSet, geneQuarSet]


v3gene = venn3_unweighted(geneSets, ('D1', 'D0.5', 'D0.25'))
plt.title('gene')
plt.savefig('genusVenn.png')
plt.clf()
plt.cla()



