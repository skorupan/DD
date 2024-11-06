import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-pt','--project_name',required=True,help='Name of project')
parser.add_argument('-fp','--file_path',required=True,help='Path to project folder without name of project folder')
parser.add_argument('-it','--n_iteration',required=True,help='Number of current iteration')

io_args = parser.parse_args()
import time

protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)

# Dictionary to store unique molecules across previous iterations
old_dict = {}
# Loop through each previous iteration up to the current one
for i in range(1,n_it):
    # Open and read training labels file from previous iteration
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/training_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1
    # Open and read validation labels file from previous iteration
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/validation_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1
    # Open and read testing labels file from previous iteration
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/testing_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1

t=time.time()
# Dictionaries to store molecules from current iteration's train, validation, and test sets
new_train = {}
new_valid = {}
new_test = {}
# Load train set molecules into `new_train` dictionary
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/train_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_train[tmpp] = 1
# Load validation set molecules into `new_valid` dictionary
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/valid_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_valid[tmpp] = 1
# Load test set molecules into `new_test` dictionary
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/test_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_test[tmpp] = 1
print(time.time()-t)

t=time.time()
# Remove duplicates between new_train, new_valid, and new_test sets
for keys in new_train.keys():
    if keys in new_valid.keys():
        new_valid.pop(keys) # Remove from validation if also in train
    if keys in new_test.keys():
        new_test.pop(keys) # Remove from test if also in train
for keys in new_valid.keys():
    if keys in new_test.keys():
        new_test.pop(keys) # Remove from test if also in validation
print(time.time()-t)

# Remove duplicates between old and new sets (seen in previous iterations)
for keys in old_dict.keys():
    if keys in new_train.keys():
        new_train.pop(keys) # Remove from train if seen in previous iterations
    if keys in new_valid.keys():
        new_valid.pop(keys) # Remove from validation if seen in previous iterations
    if keys in new_test.keys():
        new_test.pop(keys) # Remove from test if seen in previous iterations

# Write non-duplicate molecules to final train, validation, and test set files
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/train_set.txt','w') as ref:
    for keys in new_train.keys():
        ref.write(keys+'\n') # Write each molecule identifier in train set
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/valid_set.txt','w') as ref:
    for keys in new_valid.keys():
        ref.write(keys+'\n') # Write each molecule identifier in validation set
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/test_set.txt','w') as ref:
    for keys in new_test.keys():
        ref.write(keys+'\n') # Write each molecule identifier in test set
