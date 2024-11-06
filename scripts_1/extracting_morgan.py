# Reads the ids found in sampling and finds the corresponding morgan fingerprint
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name', required=True, help='Name of the DD project')
parser.add_argument('-fp', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-it', '--n_iteration', required=True, help='Number of current iteration')
parser.add_argument('-md', '--morgan_directory', required=True, help='Path to directory containing Morgan fingerprints for the database')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')

io_args = parser.parse_args()

import os
from multiprocessing import Pool
import time
from contextlib import closing
import numpy as np

protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)
morgan_directory = io_args.morgan_directory
tot_process = int(io_args.tot_process)

# Function to extract Morgan fingerprints for each molecule ID in sampling files (for model use)
def extract_morgan(file_name):
    # Initialize dictionaries to track molecules in train, valid, and test sets
    train = {}
    test = {}
    valid = {}
    # Load train, valid, and test sets for the current iteration
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/train_set.txt", 'r') as ref:
        for line in ref:
            train[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/valid_set.txt", 'r') as ref:
        for line in ref:
            valid[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/test_set.txt", 'r') as ref:
        for line in ref:
            test[line.rstrip()] = 0
    
    # Open output files to write filtered fingerprints for train, valid, and test sets
    ref1 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'train_' + file_name.split('/')[-1], 'w')
    ref2 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'valid_' + file_name.split('/')[-1], 'w')
    ref3 = open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'test_' + file_name.split('/')[-1],
                'w')
    # Read and classify each line based on molecule ID
    with open(file_name, 'r') as ref:
        flag = 0
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp in train.keys():
                train[tmpp] += 1
                fn = 1
                if train[tmpp] == 1: flag = 1
            elif tmpp in valid.keys():
                valid[tmpp] += 1
                fn = 2
                if valid[tmpp] == 1: flag = 1
            elif tmpp in test.keys():
                test[tmpp] += 1
                fn = 3
                if test[tmpp] == 1: flag = 1
            if flag == 1:
                if fn == 1:
                    ref1.write(line)
                if fn == 2:
                    ref2.write(line)
                if fn == 3:
                    ref3.write(line)
            flag = 0

# Function to read all lines from a file and return as a list
def alternate_concat(files):
    to_return = []
    with open(files, 'r') as ref:
        for line in ref:
            to_return.append(line)
    return to_return

# Function to delete a file
def delete_all(files):
    os.remove(files)

# Function to remove duplicate molecules in a Morgan file
def morgan_duplicacy(f_name):
    flag = 0
    mol_list = {} # Dictionary to track unique molecules
    ref1 = open(f_name[:-4] + '_updated.csv', 'a') # Open updated file to write unique molecules
    with open(f_name, 'r') as ref:
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp not in mol_list:
                mol_list[tmpp] = 1
                flag = 1
            if flag == 1:
                ref1.write(line) # Write unique line to updated file
                flag = 0
    os.remove(f_name)


if __name__ == '__main__':
    # Create directory to store Morgan fingerprints if it doesn't exist
    try:
        os.mkdir(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan')
    except:
        pass

    files = []
    for f in glob.glob(morgan_directory + "/*.txt"):
        files.append(f)
        
     # Parallelize extraction of fingerprints for train, valid, and test sets
    t = time.time()
    with closing(Pool(np.min([tot_process, len(files)]))) as pool:
        pool.map(extract_morgan, files)
    print(time.time() - t)

    # Process each set type (train, valid, test) for concatenation
    all_to_delete = []
    for type_to in ['train', 'valid', 'test']:
        t = time.time()
        files = []
        for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '*'):
            files.append(f)
            all_to_delete.append(f)
        print(len(files))
        if len(files) == 0:
            print("Error in address above")
            break
        with closing(Pool(np.min([tot_process, len(files)]))) as pool:
            to_print = pool.map(alternate_concat, files)
        with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '_morgan_1024.csv',
                  'w') as ref:
            for file_data in to_print:
                for line in file_data:
                    ref.write(line)
        to_print = []
        print(type_to, time.time() - t)
        
    # Get list of all concatenated files with "morgan" in filename for deduplication
    f_names = []
    for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/*morgan*'):
        f_names.append(f)
        
    # Deduplicate each concatenated file in parallel
    t = time.time()
    with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
        pool.map(morgan_duplicacy, f_names)
    print(time.time() - t)

    # Delete all temporary files generated in the process
    with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
        pool.map(delete_all, all_to_delete)
