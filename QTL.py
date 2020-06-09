def read_file():
    """
    Opens a file, puts the marker together with the bands in a dictionary.
    :return: Dictionary containing markers and corresponding the bands.
    """
    file = open("CvixLer-MarkerSubset-LG1.txt", "r")
    for _ in range(7):
        file.readline()
    alleles = {}
    allele_name = ""

    for line in file:
        if line.startswith("  "):
            bands = line.strip().split(" ")
            for band in bands:
                if band != "":
                    alleles[allele_name].append(band)
        else:
            allele_name = line.strip().split("(")[0]
            alleles[allele_name] = []
    return alleles


def compare_alleles(alleles):
    """
    Compares the bands of different markers and counts the number of recombinants
    and the total number of compared bands.
    :param alleles: Dictionary containing markers and corresponding the bands
    :return: A dictionary with keys containing compared markers.
    The value contains a list with the number of recombinants and the total number of compared bands.
    """
    compared = {}
    for allele1, val1 in alleles.items():
        for allele2, val2 in alleles.items():
            if allele1 != allele2:
                compared[allele1 + "\t" + allele2] = compare(val1, val2)
    return compared


def compare(allele1, allele2):
    """
    Compares two alleles to find recombinants.
    :param allele1: String containing band.
    :param allele2: String containing band.
    :return: A list containing the number of recombinants and the total number of compared bands.
    """
    total = 0
    recombinants = 0
    for i in range(len(allele1)):
        if allele1 != "-" and allele2 != "-":
            total += 1
            if allele1[i] != allele2[i]:
                recombinants += 1
    return [recombinants, total]


def calculate_distance(compared_alleles):
    """
    Firstly determines the highest found distance of each marker.
    Secondly subtracts the lowest number from the highest found number from each of the highest found numbers.
    Thridly sorts the numbers from low to high.
    :param compared_alleles: A dictionary with keys containing compared markers.
    The value contains a list with the number of recombinants and the total number of compared bands.
    :return: Dictionary with marker names as keys and marker distances as values.
    """
    allele_distance = {}
    allele = None
    highest = int
    lowest = 100
    for key, value in compared_alleles.items():
        distance = round(value[0] / value[1] * 100, 2)
        if allele != key.split("\t")[0]:
            if allele is not None:
                if highest < lowest:
                    lowest = highest
                allele_distance[allele] = highest
            allele = key.split("\t")[0]
            highest = distance
        else:
            if distance > highest:
                highest = distance

    if highest < lowest:
        lowest = highest
    allele_distance[allele] = highest

    return adjust_alleles(allele_distance, lowest)


def adjust_alleles(allele_distance, lowest):
    """
    Subtracts the lowest number from the highest found number from each of the highest found numbers
    and sorts them from low to high.
    :param allele_distance: Dictionary with marker names as keys and marker distances as values.
    :param lowest: int containing the lowest of the highest marker distances found.
    :return: Dictionary with marker names as keys and marker distances as values.
    """
    allele_dist_ajdust = {}
    for key, value in allele_distance.items():
        allele_dist_ajdust[key] = value - lowest

    sorted_alleles = sorted(allele_dist_ajdust.items(), key=lambda x: x[1])
    return sorted_alleles


def save_to_file(sorted_alleles):
    """
    Saves markers with distances to a file.
    :param sorted_alleles: Dictionary with marker names as keys and marker distances as values.
    """
    file = open("QTL output.txt", "w")
    for allele, cm in sorted_alleles:
        file.write(allele + "\t" + str(round(cm, 2)) + "\n")
    file.close()


if __name__ == '__main__':
    found_alleles = read_file()
    alleles_compared = compare_alleles(found_alleles)
    alleles_sorted = calculate_distance(alleles_compared)
    save_to_file(alleles_sorted)
