# Docking scores for each molecular set are saved to a corresponding txt file inside iteration_1. The number of processes should match the number of docking sdf files (three). The score_keyword argument must match the title of the field that stores the docking score of a molecule in the sdf files.
import argparse
import builtins as __builtin__
import glob
import gzip
import os
from contextlib import closing
from multiprocessing import Pool


# For debugging purposes only:
def print(*args, **kwargs):
    __builtin__.print('\t extract_L: ', end="")
    return __builtin__.print(*args, **kwargs)


parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name', required=True, help='Name of the DD project')
parser.add_argument('-fp', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-n_it', '--iteration_no', required=True, help='Number of current iteration')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')
parser.add_argument('-score', '--score_keyword', required=True, help='Score keyword. Name of the field storing the docking score in the SDF files of docking results')

io_args = parser.parse_args()
n_it = int(io_args.iteration_no)
protein = io_args.project_name
file_path = io_args.file_path
tot_process = int(io_args.tot_process)
key_word = str(io_args.score_keyword)

# mol_key = 'ZINC'
print("Keyword: ", key_word)

# Function to extract docking scores (labels) from each molecule in the file
def get_scores(ref):
    scores = []
    for line in ref:  # Looping through the molecules
        zinc_id = line.rstrip() # Store the molecule's unique ID. *Note: zinc_id represents molecule ID
        line = ref.readline()
        # '$$$' signifies end of molecule info
        while line != '' and line[:4] != '$$$$':  # Looping through its information and saving scores

            tmp = line.rstrip().split('<')[-1] # Extract the keyword field name

            if key_word == tmp[:-1]:
                tmpp = float(ref.readline().rstrip())
                # Check if the score is within a reasonable range, otherwise print it for debugging
                if tmpp > 50 or tmpp < -50:
                if tmpp > 50 or tmpp < -50:
                    print(zinc_id, tmpp) # Append the score and ID to the list
                else:
                    scores.append([zinc_id, tmpp])

            line = ref.readline()
    return scores

# Function to extract scores from a compressed docking file
def extract_glide_score(filen):
    scores = []
    try:
        # Opening the GNU compressed file
        with gzip.open(filen, 'rt') as ref:
            scores = get_scores(ref)

    except Exception as e:
        print('Exception occured in Extract_labels.py: ', e)
        # file is already decompressed
        with open(filen, 'r') as ref:
            scores = get_scores(ref)

    # Determine the dataset type (training, testing, or validation) based on the filename
    if 'test' in os.path.basename(filen):
        new_name = 'testing'
    elif 'valid' in os.path.basename(filen):
        new_name = 'validation'
    elif 'train' in os.path.basename(filen):
        new_name = 'training'
    else:
        print("Could not generate new training set")
        exit()
        
    # Write the extracted scores to a labeled output file
    with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/' + new_name + '_' + 'labels.txt', 'w') as ref:
        ref.write('r_i_docking_score' + ',' + 'ZINC_ID' + '\n')
        for z_id, gc in scores:
            ref.write(str(gc) + ',' + z_id + '\n')


if __name__ == '__main__':
    files = []  # List to hold paths of all docking files
    iter_path = file_path + '/' + protein + '/iteration_' + str(n_it)

    # Checking to see if the labels have already been extracted:
    sets = ["training", "testing", "validation"]
    files_labels = glob.glob(iter_path + "/*_labels.txt")
    foundAll = True
    for s in sets:
        found = False
        print(s)
        for f in files_labels:
            set_name = f.split('/')[-1].split("_labels.txt")[0]
            if set_name == s:
                found = True
                print('Found')
                break
        if not found:
            foundAll = False
            print('Not Found')
            break
    if foundAll:
        print('Labels have already been extracted...')
        print('Remove "_labels.text" files in \"' + iter_path + '\" to re-extract')
        exit(0)
        
    # Glob pattern to find all docking result files in the current iteration directory
    path = iter_path + '/docked/*.sdf*'
    path_labels = iter_path + '/*labels*'

    # Append each file found to the list of files
    for f in glob.glob(path):
        files.append(f)

    print("num files in", path, ":", len(files))
    print("Files:", [os.path.basename(f) for f in files])
    if len(files) == 0:
        print('NO FILES IN: ', path)
        print('CANCEL JOB...')
        exit(1)

    # Parallel running of the extract_glide_score() with each file path of the files array
    with closing(Pool(len(files))) as pool:
        pool.map(extract_glide_score, files)
