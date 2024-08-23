import os


def list_input_creation(filedirectory, outputfile):
    os.chdir(filedirectory)

    with open(outputfile, 'a') as handle:
        for file in os.listdir():
            if file.endswith('.bam.fasta'):
                handle.write(str(file) + '\n')
    handle.close()


list_input_creation('/home/albertotr/OneDrive/Data/'
                    'Cambridge_Project/Mapped_output_VA94_7994_1_7P/', 'List_files_VA94_7994_1_7P.txt')
list_input_creation('/home/albertotr/OneDrive/Data/'
                    'Cambridge_Project/Mapped_output_Rlow/', 'List_files_Rlow.txt')
