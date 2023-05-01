import os
import install_packages
import argparse


def run_blast(indir, outdir, reference):
    for genome in os.listdir(indir):
        subject_file = os.path.join(indir, genome)
        blast_out = os.path.join(outdir, f"{genome[:-6]}_tmp.fasta")

        print(f"Running blastn for", f"query {reference}", f"subject {subject_file}.", sep='\n')
        os.system(f"blastn -query {reference} -subject {subject_file} -outfmt \"6 delim=@ sseq\" -out {blast_out}")

        with open(subject_file) as get_name:
            name = get_name.readline()
        with open(blast_out, 'r') as sequence:
            out_file = open(os.path.join(outdir, f"{genome[:-6]}.fasta"), 'w')
            out_file.write(name)
            out_file.write(sequence.read())
            out_file.close()
        os.system(f"rm {blast_out}")


def join_genes(dir):
    selfname = "general.fasta"
    with open(os.path.join(dir, selfname), 'w') as megafasta:
        for fasta_file in os.listdir(dir):
            if fasta_file != selfname:
                with open(os.path.join(dir, fasta_file), 'r') as fasta:
                    megafasta.write(fasta.read())


def run_mafft(outdir, infile):
    out_name = "mafft_results.fasta"
    os.system(f"mafft {infile} > {os.path.join(outdir, out_name)}")


def main(indir, outdir, reference):
    install_packages.main()

    os.system(f"mkdir {outdir}")

    run_blast(indir, outdir, reference)
    join_genes(outdir)
    run_mafft(outdir, os.path.join(outdir, "general.fasta"))

    print("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="geneMSA",
        description="Multiple sequence alignment for proteins from given genomes and reference.\n"
                    "Note: input and reference files must be FASTA!")

    parser.add_argument('-i', "--input-files-directory")
    parser.add_argument('-o', "--output-files-directory", default="geneMSA_output")
    parser.add_argument('-r', "--reference")

    args = vars(parser.parse_args())
    input_files_directory = args["input_files_directory"]
    output_files_directory = args["output_files_directory"]
    reference = args["reference"]
    print("Multiple sequence alignment for proteins from given genomes and reference.",
          f"-i {input_files_directory}",
          f"-o {output_files_directory}",
          f"-r {reference}", sep='\n')

    main(input_files_directory, output_files_directory, reference)
