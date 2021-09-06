#!/usr/bin/env python

'''

    Author: Sofia Papadimitriou

    == INPUT == : the "1KGP_training_features.txt" file with the 1KGP training sets (200 combinations per set) and the
    "DIDA_training_features.txt" file with the DIDA training set (200 combinations)


    == PROCESS == : this script takes the input training data, creates 500 balanced sets and trains
    500 Random Forest predictors


    == OUTPUT == : a cPickle (.p) file with the trained model in the form of a LIST of 500 Random Forests
                   The output file is ready to be used in VarCoPP.py.


     Usage examples
   ===========================================================================
   You can directly run the script with all default parameters, i.e. the ones described in the paper):
   $ python VarCoPP_train.py

   In case you want to use different training set files, make sure to provide their complete path and name.

   In case you want to use a different prefix for the output file, you can specify it in the option -f.
   $ python VarCoPP_train.py -f test

   Then, the output model will be named as: test_model.p

   You can also specify a different output directory, which will be created in case it doesn't exist:
   $ python VarCoPP_train.py -o test_directory


   '''

from __future__ import division
from sklearn.ensemble import RandomForestClassifier
import os.path
import argparse
import random
import pickle


def input_parser():
    """ Function to construct and parse input parameters """

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-n', '--neutral_dataset', help='.txt file with the 1KGP training features (DEFAULT: ./training_data/1KGP_training_features.txt)',
                        default = './training_data/1KGP_training_features.txt')
    parser.add_argument('-d', '--dida_dataset',
                        help='.txt file with the DIDA training features (DEFAULT: ./training_data/DIDA_training_features.txt',
                        default = './training_data/DIDA_training_features.txt')
    parser.add_argument('-o', '--outdir', help='Output directory - full path (DEFAULT: ./)', default = './')
    parser.add_argument('-f', '--filename',
                        help='Specify the prefix of the file before the _symmary.txt and ' +
                             '_iterations.txt (DEFAULT: training)', default='training')
    parser.add_argument('-e', '--features',
                        help='list of the feature names used, delimited by comma (do NOT provide it, if you wish to re-train the' +
                             ' data with the features used in the paper)',
                        default='Flex1,Hydr1,CADD1,CADD2,CADD3,CADD4,HI_A,RecA,HI_B,RecB,Biol_Dist')
    parser.add_argument('-t', '--num_trees',
                        help='number of Decision Trees for each random forest (DEFAULT: 100)',
                        default='100')
    parser.add_argument('-m', '--max_depth',
                        help='maximum depth of each Decision Tree (DEFAULT: 10)',
                        default='10')
    return vars(parser.parse_args())


def feat_to_vector(al, feat_list):

    '''
    Function that takes a list of features and transforms them into computer
    readable values'
    '''

    # initiate the vector
    vector = []
    h = 0  # vector index

    # initiate dictionary with the vector start position of each feature
    feat_pos_dict = {}

    # iterate over the features to vectorize them
    for i in range(len(feat_list)):

        # vectorize for the aminoacid property change
        if feat_list[i] in ['Flex1', 'Hydr1']:

            # Categorize value or wt
            try:
                if al[i] in ['wt', 'wild_type']:
                    vector += [0.0]
                elif al[i] == 'NA':
                    vector += [vector[-1]]
                else:
                    vector += [float(al[i])]
            except ValueError:
                print('##### Some values of Hydr1 or Flex1 are invalid.' + \
                ' Please use the calc_aa_diff.py to annotate variants' + \
                ' with Hydr and Flex amino acid differences.')
                print('##### Script is terminated')

                quit()

            # save vector position of the feature
            feat_name = feat_list[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

        # vectorize for CADD
        elif feat_list[i] in ['CADD1', 'CADD2', 'CADD3', 'CADD4']:

            # Categorize value or wt
            try:
                if al[i] in ['wt', 'wild_type']:
                    vector += [-3.0]
                elif al[i] == 'NA':
                    vector += [vector[-1]]
                else:
                    vector += [float(al[i])]
            except ValueError:
                print('##### Some CADD raw score values are invalid.' + \
                ' Please provide a score for CADD for all variants.')
                print('##### Script is terminated')

                quit()

            # save vector position of the feature
            feat_name = feat_list[i]
            feat_pos_dict[feat_name] = [h, len(vector) - 1]
            h = len(vector)

    ############################ collect information on genes and pair info

    # iterate over the gene and pair features
    for i in range(len(feat_list)):

        ############# Haploinsufficiency probability
        if feat_list[i] in ['HI_A', 'HI_B']:

            try:
                if al[i] in ['N/A', 'NA', 'NaN', 'nan']:
                    vector += [0.19898]
                else:
                    vector += [float(al[i])]
            except ValueError:
                print('##### Some gene haploinsufficiency values are invalid.' + \
                ' Please provide a haploinsufficiency value or write' + \
                ' "NA","NaN", or "nan".')
                print('##### Script is terminated')

                quit()

            # save vector position of the feature
            feat_pos_dict[feat_list[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        ############# Recessiveness probability
        elif feat_list[i] in ['RecA', 'RecB']:

            try:
                if al[i] in ['N/A', 'NA', 'NaN', 'nan']:
                    vector += [0.12788]
                else:
                    vector += [float(al[i])]
            except ValueError:
                print('##### Some gene recessiveness values are invalid.' + \
                ' Please provide a recessiveness value or write' + \
                ' "NA","NaN", or "nan".')
                print('##### Script is terminated')

                quit()

            # save vector position of the feature
            feat_pos_dict[feat_list[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

        ############ Biological distance
        elif feat_list[i] in ['Biol_Dist']:

            try:
                if al[i] in ['N/A', 'NA', 'NaN', 'nan']:
                    vector += [11.44]
                else:
                    vector += [float(al[i])]
            except ValueError:
                print('##### Some Biological distance values are invalid.' + \
                ' Please provide a value or write' + \
                ' "NA","NaN", or "nan".')
                print('##### Script is terminated')

                quit()

            # save vector position of the feature
            feat_pos_dict[feat_list[i]] = [h, len(vector) - 1]
            # initiate vector index
            h = len(vector)

    return vector, feat_pos_dict



if __name__ == '__main__':

    # Calling the parser function to start parsing the argument
    parser = input_parser()

    ############################## check existence of output folder and correct

    output_dir = parser['outdir']
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

        print('##### Output directory did not exist, it was created during the process')

    if not output_dir.endswith('/'):
        output_dir += '/'


    ######################################### get features
    print('##### Getting features')
    feat_list = parser['features'].split(',')


    ######################################### collect DIDA features set and transform data
    print('##### Collecting the DIDA features data set')

    # collect DIDA training set
    try:
        dida_setfile = open(parser['dida_dataset'])
        dida_setfile.readline()
    except IOError:
        print('##### The DIDA training features set could not be uploaded. Make sure that' + \
        ' the training file is in the same directory with the script. ')
        print('##### Script is terminated')
        quit()

    # iterate over the combinations and convert combinations into vectors
    dida_set = {}

    for line in dida_setfile:
        line_list = line.rstrip('\n').split('\t')
        comb_id = line_list[0]

        # convert to vectors
        initial_vector = line_list[1:]  # remove the ID column
        vector, feat_order = feat_to_vector(initial_vector, feat_list)

        #add information
        dida_set[comb_id] = {'comb': vector, 'class':1}

    dida_setfile.close()


    ######################################### collect 1KGP features set and transform data
    print('##### Collecting the 1KGP features data set')

    # collect 1KGP training set
    try:
        neutral_setfile = open(parser['neutral_dataset'])
        neutral_setfile.readline()
    except IOError:
        print('##### The 1KGP training features set could not be uploaded. Make sure that' + \
        ' the training file is in the same directory with the script. ')
        print('##### Script is terminated')
        quit()

    # iterate over the combinations and convert combinations into vectors
    neutral_set = []
    set_combs = {}  # initialize the dictionary to hold combinations of a set

    set_ids = [] #collect all set IDs

    for line in neutral_setfile:
        line_list = line.rstrip('\n').split('\t')

        set_id = line_list[0].split('_')[:2] #collect the ID of the set
        comb_id = line_list[0] #collect the ID of the combination

        #if we reach a new set, add information and save the combinations of the previous one
        if set_id not in set_ids:

            #save the information of the previous set
            if set_combs!={}:
                neutral_set += [set_combs]

            #add the ID in the list
            set_ids += [set_id]

            # convert to vectors and re-initialize
            initial_vector = line_list[1:]  # remove the ID column
            vector, feat_order = feat_to_vector(initial_vector, feat_list)
            set_combs = {comb_id : {'comb': vector, 'class':0}}

        #otherwise, add information to the current set
        else:
            # convert to vectors
            initial_vector = line_list[1:]  # remove the ID column
            vector, feat_order = feat_to_vector(initial_vector, feat_list)
            set_combs[comb_id] = {'comb': vector, 'class':0}


    #collect the last set
    neutral_set += [set_combs]

    neutral_setfile.close()


    ######################################### train VarCoPP

    rf_models = []

    for subset in range(0, len(neutral_set)):

        print('##### Training VarCoPP with training set %d'%(subset+1))

        # define the current neutral set
        neutral_subset = neutral_set[subset]

        #create the balanced set
        balanced_set = neutral_subset.copy()
        balanced_set.update(dida_set)

        #randomize the balanced set
        balanced_set_keys = list(balanced_set.keys())
        random.shuffle(balanced_set_keys)

        #create instance and label sets
        instance_set = []
        label_set = []

        # append dictionary information in the lists
        for key in balanced_set_keys:
            instance_set += [balanced_set[key]['comb']]
            label_set += [balanced_set[key]['class']]

        # train Random Forest
        clf = RandomForestClassifier(n_estimators=int(parser['num_trees']), max_depth=int(parser['max_depth']))
        new_clf = clf.fit(instance_set, label_set)

        # add random forests for each set
        rf_models += [[new_clf]]


    ######################################### save VarCoPP model
    print('##### Saving VarCoPP model')
    with open(output_dir + parser['filename'] + '_model.p', 'wb') as random_f_file:
        pickle.dump(rf_models, random_f_file)



























