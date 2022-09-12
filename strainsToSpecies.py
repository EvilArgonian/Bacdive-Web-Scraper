inFile = "mergeTheseStrains.tsv"
# inFile = "testMergeInput.tsv"
outFile = "mergedSpecies.txt"
thresholdDecimal = 0.75

clearOutFile = open(outFile, "a")
clearOutFile.truncate(0)
clearOutFile.close()

mergeDict = {}
with open(inFile, "r") as input:
    for line in input.readlines():
        # [strain_id, species, strain_des, culture_colNo, taxDomain, taxPhylum, taxClass, taxOrder, taxFamily, taxGenus, living_state, pathogenicity, aerobicity, gc_content, tags1, tags2, iso_source] = line.split("\t")
        speciesName = line.split("\t")[1]
        strainData = line.split("\t")
        if speciesName in mergeDict.keys():
            newMergeGroup = mergeDict.get(speciesName)
            newMergeGroup.append(strainData)
            mergeDict.update({speciesName:newMergeGroup})
        else:
            mergeDict.update({speciesName:[strainData]})


for mergeKey in mergeDict:
    toBeMerged = mergeDict.get(mergeKey)
    species = toBeMerged[-1][1]
    taxDomain = toBeMerged[-1][4]
    taxPhylum = toBeMerged[-1][5]
    taxClass = toBeMerged[-1][6]
    taxOrder = toBeMerged[-1][7]
    taxFamily = toBeMerged[-1][8]
    taxGenus = toBeMerged[-1][9]

    for item in toBeMerged:
        if species != item[1] \
                or (taxDomain and item[4] and taxDomain != item[4]) \
                or (taxPhylum and item[5] and taxPhylum != item[5]) \
                or (taxClass and item[6] and taxClass != item[6]) \
                or (taxOrder and item[7] and taxOrder != item[7]) \
                or (taxFamily and item[8] and taxFamily != item[8]) \
                or (taxGenus and item[9] and taxGenus != item[9]):  # If any non-null entry doesn't match
            print("Tried to merge a group including a differing species! Failing!")
            print(
                species + "\t" + taxDomain + "\t" + taxPhylum + "\t" + taxClass + "\t" + taxOrder + "\t" + taxFamily + "\t" + taxGenus)
            print(item[1] + "\t" + item[4] + "\t" + item[5] + "\t" +
                  item[6] + "\t" + item[7] + "\t" + item[8] + "\t" +
                  item[9])
            exit(1)

    pathogenicity = "NULL"
    for item in toBeMerged:
        thisPathogenicity = item[11]
        if thisPathogenicity == "True":  # Supercedes any others
            pathogenicity = "Pathogenic"
            break
        if thisPathogenicity == "Infection" and pathogenicity != "Pathogenic":  # Supercedes Not Pathogenic
            pathogenicity = "Infection"
        if thisPathogenicity == "False" and pathogenicity == "NULL":  # Only occurs if still NULL
            pathogenicity = "Not Pathogenic"

    living_state = "NULL"
    if pathogenicity == "Pathogenic":
        living_state = "Not Free Living"
    else:
        counts = {}
        for item in toBeMerged:
            thisLivingState = item[10]
            counts.update({"NULL": 0})  # Guarantees function of discounting NULLs
            counts.update({"": 0})  # Also discounting blanks, which seems to also be possible
            if thisLivingState in counts.keys():
                counts.update({thisLivingState: counts.get(thisLivingState) + 1})
            else:
                counts.update({thisLivingState: 1})
        totalNonNullEntries = sum(counts.values()) - counts.get("NULL") - counts.get("")
        threshold = totalNonNullEntries * thresholdDecimal
        for key in counts.keys():
            if counts.get(key) > threshold:  # Supercedes any others if clear threshold
                living_state = key
                break
            if key == "Not Free Living":  # Supercedes NULL and Free Living
                living_state = "Not Free Living"
            if living_state == "NULL" and key == "Free Living":  # Only occurs if still NULL
                living_state = "Free Living"

    aerobicity = "NULL"
    counts = {}
    for item in toBeMerged:
        thisAerobicity = item[12].lower().strip()
        # print("Aerobicity: " + thisAerobicity + "\n")
        counts.update({"null": 0})  # Guarantees function of discounting NULLs
        counts.update({"": 0})  # Also discounting blanks, which seems to also be possible
        if thisAerobicity in counts.keys():
            counts.update({thisAerobicity: counts.get(thisAerobicity) + 1})
        else:
            counts.update({thisAerobicity: 1})
    totalNonNullEntries = sum(counts.values()) - counts.get("null") - counts.get("")
    threshold = totalNonNullEntries * thresholdDecimal
    for key in counts.keys():
        if counts.get(key) > threshold:  # Supercedes any others if clear threshold
            aerobicity = key
            break
        if key == "anaerobic":  # Supercedes any non-threshold-breakers
            aerobicity = "anaerobic"
        if aerobicity != "anaerobic" and key == "obligate anaerobe":  # Supercedes NULL, aerobe, obligate aerobe, microaerophile, and facultative anaerobe
            aerobicity = "obligate anaerobe"
        if aerobicity != "anaerobic" and aerobicity != "obligate anaerobe" and key == "microaerophile":  # Supercedes NULL, aerobe, obligate aerobe, and facultative anaerobe
            aerobicity = "microaerophile"
        if aerobicity != "anaerobic" and aerobicity != "obligate anaerobe" and aerobicity != "microaerophile" and key == "facultative anaerobe":  # Supercedes NULL, aerobe, and obligate aerobe
            aerobicity = "facultative anaerobe"
        if aerobicity != "anaerobic" and aerobicity != "obligate anaerobe" and aerobicity != "microaerophile" and aerobicity != "facultative anaerobe" and key == "obligate aerobe":  # Supercedes NULL and aerobe
            aerobicity = "obligate aerobe"
        if aerobicity == "NULL" and key == "aerobe":  # Only occurs if still NULL
            aerobicity = "aerobe"

    gc_content = 0.0
    numApproxGCs = 0
    sumApproxGCs = 0.0
    for item in toBeMerged:
        if item[13] != "NULL" and item[13] != "":
            numApproxGCs += 1
            sumApproxGCs += float(item[13])
    if numApproxGCs > 0:
        gc_content = float(
            sumApproxGCs / numApproxGCs)  # Overwrites with average of strain GC contents, if extant

    with open(outFile, "a+") as out:
        out.write(species + "\t" +
                  taxDomain + "\t" +
                  taxPhylum + "\t" +
                  taxClass + "\t" +
                  taxOrder + "\t" +
                  taxFamily + "\t" +
                  taxGenus + "\t" +
                  living_state + "\t" +
                  pathogenicity + "\t" +
                  aerobicity + "\t" +
                  "{:.3f}".format(gc_content) + "\n")

    print(species + "\t" +
          taxDomain + "\t" +
          taxPhylum + "\t" +
          taxClass + "\t" +
          taxOrder + "\t" +
          taxFamily + "\t" +
          taxGenus + "\t" +
          living_state + "\t" +
          pathogenicity + "\t" +
          aerobicity + "\t" +
          "{:.3f}".format(gc_content))
