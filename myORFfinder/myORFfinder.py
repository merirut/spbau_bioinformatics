import argparse


def detect_amino(kmer):
    aminos = {
        "AGG": "R", "AGA": "R", "AGC": "S", "AGT": "S",  # arginine, serine
        "AAG": "K", "AAA": "K", "AAC": "N", "AAT": "N",  # lysine, asparagine
        "ACG": "T", "ACA": "T", "ACC": "T", "ACT": "T",  # threonine
        "ATG": "M", "ATA": "I", "ATC": "I", "ATT": "I",  # methionine (start), isoleucine

        "CGG": "R", "CGA": "R", "CGC": "R", "CGT": "R",  # arginine
        "CAG": "Q", "CAA": "Q", "CAC": "H", "CAT": "H",  # glutamine, histidine
        "CCG": "P", "CCA": "P", "CCC": "P", "CCT": "P",  # proline
        "CTG": "L", "CTA": "L", "CTC": "L", "CTT": "L",  # leucine

        "TGG": "W", "TGA": "stop", "TGC": "C", "TGT": "C",  # tryptophan, stop, cysteine
        "TAG": "stop", "TAA": "stop", "TAC": "Y", "TAT": "Y",  # stop, stop, tyrosine
        "TCG": "S", "TCA": "S", "TCC": "S", "TCT": "S",  # serine
        "TTG": "L", "TTA": "L", "TTC": "F", "TTT": "F",  # leucine, phenylalanine

        "GGG": "G", "GGA": "G", "GGC": "G", "GGT": "G",  # glycine
        "GAG": "E", "GAA": "E", "GAC": "D", "GAT": "D",  # glutamate, aspartate
        "GCG": "A", "GCA": "A", "GCC": "A", "GCT": "A",  # alanine
        "GTG": "V", "GTA": "V", "GTC": "V", "GTT": "V"  # valine
    }

    return aminos[kmer]


def get_sequence(file_path):
    mRNAfile = open(file_path, 'r')
    raw = mRNAfile.readlines()
    mRNAfile.close()

    paired_list = []
    name = None
    seq = ""

    for line in raw:
        if (line[0] == '>'):
            if name:
                paired_list.append([name, seq])
            name = line[:-1]
            seq = ""
        else:
            seq += line[:-1]
    paired_list.append([name, seq])
    return paired_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="myORFfinder",
        description="Finds all Open Reading Frames in a given file with sequences.")

    parser.add_argument("filename")
    parser.add_argument('-o', "--output", type=str, default="myORFfinder_output.txt", help="Output name")
    parser.add_argument('-l', "--min-length", type=int, default=30, help="Min ODF length (default: 30)")

    args = vars(parser.parse_args())

    pairs = get_sequence(args["filename"])
    outfile = open(args["output"], 'w')
    min_length = args["min_length"]

    for label, mRNA in pairs:

        outfile.write(label + "\n")

        mRNAlen = len(mRNA)
        start = mRNA.find("ATG") # first start-codon

        while (start < mRNAlen - 2):
            CDS = ""
            aminoseq = ""
            CDS_started = False
            CDS_finished = False
            current_start = start

            for i in range(current_start, mRNAlen - 2, 3):
                kmer = mRNA[i:i + 3]

                if (kmer == "ATG"): # start-codon
                    CDS_started = True
                    start = i + 3

                if (kmer == "TAA" or kmer == "TAG" or kmer == "TGA"): # stop-codon
                    CDS_finished = True
                    break

                if (CDS_started and not CDS_finished):
                    CDS += kmer
                    aminoseq += detect_amino(kmer)


            if (CDS_finished and len(CDS) >= min_length):
                outfile.write(CDS + "\n" + aminoseq + "\n")

            if (not CDS_started):
                start += 1

    outfile.close()
