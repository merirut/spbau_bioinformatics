import os


def install_blast():
    os.system("conda create -n blast && conda activate blast")
    os.system("conda install -c bioconda blast==2.12.0 --yes")


def install_mafft():
    os.system("conda install -c bioconda --yes mafft")


def main():
    if len(os.popen("conda list blast").readlines()) < 3:
        install_blast()
    #if len(os.popen("conda list mafft").readlines()) < 3:
        #install_mafft()


if __name__ == "__main__":
    main()
