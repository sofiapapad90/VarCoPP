#!/usr/bin/env python

'''

   Author: Sofia Papadimitriou

   == INPUT == : the "1KGP_training_features.txt" file with the 1KGP training sets (200 combinations per set),
     the "DIDA_training_features.txt" file with the DIDA training set (200 combinations) and the "training_data_pairs.txt"
     file with the pair information for all combination IDs present in the training data.


    == PROCESS == : this script takes the input training data, creates 500 balanced sets and creates a
    cross-validation procedure in each one of them to evaluate the model's performance


    == CROSS-VALIDATION == : it is a LeaveOnePairOut method that takes place inside
    EACH individual balanced set


    == OUTPUT == : 1. 'training_cross_val_summary.txt', with a summary of the average performance among all 500
                   individual predictors
                   2. 'training_cross_val_iterations.txt', with a detailed performance for each of the 500
                   individual predictors
                   3. 'training_metrics_dict.json', with statistics that can be used later for plots
                   4. 'feat_dict.json', with the importance of all features between all 500 individual predictors


    Usage examples
   ===========================================================================
   You can directly run the script with all the default parameters (the ones used for the paper as well):
   Usage: $ python VarCoPP_cross_validation.py

   In case you want to use a different output filename, you can specify it in the option -f.
   Usage: $ python VarCoPP_cross_validation.py -f my_filename

   In case you want to change the disease-causing probability threshold to the standard 0.50, you can use the
   -s CLASS_THRESHOLD parameter.
   Usage: $ python VarCoPP_cross_validation.py -s 0.5

   *** NOTE:  You don't have to specify the -e FEATURES parameter, unless you want to remove some features from the
   features list.


'''

from __future__ import division
from sklearn.ensemble import RandomForestClassifier  # for Random Forest functions
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn import metrics  # for statistics
from sklearn.model_selection import LeaveOneGroupOut
import math
import numpy as np
import os.path
import argparse
import json
import re



def input_parser():
    """ Function to construct and parse input parameters """

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-n', '--neutral_dataset', help='.txt file with the 1KGP training features (DEFAULT: 1KGP_training_features.txt)',
                        default='./1KGP_training_features.txt')
    parser.add_argument('-d', '--dida_dataset',
                        help='.txt file with the DIDA training features (DEFAULT: DIDA_training_features.txt)',
                        default='./DIDA_training_features.txt')
    parser.add_argument('-p', '--training_pairs',
                        help='.txt file with gene pair information of all training data (DEFAULT: training_data_pairs.txt',
                        default = './training_data_pairs.txt')
    parser.add_argument('-o', '--outdir', help='Output directory - full path (DEFAULT: ./)', default = './')
    parser.add_argument('-f', '--filename',
                        help='Specify the name of the file before the _symmary.txt and ' +
                             '_iterations.txt (DEFAULT: training_cross_val)', default='training_cross_val')
    parser.add_argument('-e', '--features',
                        help='list of the feature names used, delimited by comma (do NOT provide it if you wish to re-train the' +
                             ' data with the features used in the paper)',
                        default='Flex1,Hydr1,CADD1,CADD2,CADD3,CADD4,HI_A,RecA,HI_B,RecB,Biol_Dist')
    parser.add_argument('-t', '--num_trees',
                        help='number of Decision Trees for each random forest (DEFAULT: 100)',
                        default='100')
    parser.add_argument('-m', '--max_depth',
                        help='maximum depth of each Decision Tree (DEFAULT: 10)',
                        default='10')
    parser.add_argument('-s', '--class_threshold',
                        help='Specify a minimum probability threshold for class 1 (DEFAULT: 0.489)',
                        default='0.489')

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


def mean(the_list):
    '''Function that takes a list with numerical elements and returns the mean
    value of the elements'''

    total = 0
    for i in the_list:
        total += i
    mean = total / len(the_list)

    return mean


def sd(the_list, mean):
    '''Function that takes a list with numerical elements and returns the standard
    deviation of the element values'''

    sum_of_sq = 0

    # calculate the sum of (value-mean)**2
    for value in the_list:
        sum_of_sq += (value - mean) ** 2

        # calculate variance
    variance = sum_of_sq / (len(the_list))

    # calculate standard deviation

    sd = math.sqrt(variance)

    return sd


def sorted_strings(l):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def stratified_gene_pair_groups(set_dict, instance_keys):
    ''' Function that takes a dictionary (set_dict) with keys and all information about
    instances (including the key 'pair' that contains a tuple with the pair of genes),
    and a list of the keys that correspond to instances (in the same order as the instances)
    are present in the set.

    Returns an array that corresponds to numbers based on the unique gene pair each instance
    belongs to, in the same order as these instances appear in the instance set'''

    # list with pair group number
    paired_groups = []

    # pair_group number dictionary
    pair_num_dict = {}

    # iterate over the instances
    for i in range(len(instance_keys)):
        key = instance_keys[i]
        # add pairs in dictionaries
        geneA, geneB = set_dict[key]['pair']
        pair_num_dict[(geneA, geneB)] = 0

    ############ enumerate pairs
    i = 1
    for key in pair_num_dict:
        pair_num_dict[key] = i

        i += 1

    # create dictionary with keys and group number
    for key in instance_keys:
        geneA, geneB = set_dict[key]['pair']
        paired_groups += [pair_num_dict[(geneA, geneB)]]

    return paired_groups


def predict_classes(class_probs, threshold):
    '''Function that takes the probabilities of class assignments for each
    element and the threshold for class 1 and returns an array with class
    labels for each element in the list'''

    # list to hold labels
    predicted_classes = []

    # iterate over the elements
    for el in class_probs:
        # check threshold for class 1
        if el[1] > threshold:
            predicted_classes += [1]
        else:
            predicted_classes += [0]

    predicted_classes = np.array(predicted_classes)

    return predicted_classes


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


    ######################################### collect pairs
    print('##### Collecting training data pairs (for cross-validation)')
    training_pairs = {}

    pairsfile = open(parser['training_pairs'])
    pairsfile.readline()

    for line in pairsfile:
        line_list = line.rstrip('\n').split('\t')

        training_pairs[line_list[0]] = [line_list[1], line_list[2]]

    pairsfile.close()


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
        dida_set[comb_id] = {'comb': vector, 'class':1, 'pair':training_pairs[comb_id]}

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

            set_combs = {comb_id : {'comb': vector, 'class':0, 'pair':training_pairs[comb_id]}}

        #otherwise, add information to the current set
        else:
            # convert to vectors
            initial_vector = line_list[1:]  # remove the ID column
            vector, feat_order = feat_to_vector(initial_vector, feat_list)
            set_combs[comb_id] = {'comb': vector, 'class':0, 'pair':training_pairs[comb_id]}



    #collect the last set
    neutral_set += [set_combs]

    neutral_setfile.close()


    ####################################### Collect balanced sets in one big file
    print('##### Creating balanced training sets')
    set_dicts = []
    for subset in neutral_set:

        subset_dict = subset.copy()
        subset_dict.update(dida_set)

        set_dicts += [subset_dict]


    ######################################################### open output files
    iterations_out_file = open(output_dir + parser['filename'] + '_iterations.txt', 'w')
    summary_out_file = open(output_dir + parser['filename'] + '_summary.txt', 'w')


    ###################################### initiate combination number and vector length for statistics
    print('##### Preparing statistics')

    combs_num = len(list(set_dicts[0].keys()))
    vector_length = len(feat_list)

    ##################################### feature importance for statistics

    feat_imp_dict = {}
    for feat in feat_order:
        feat_imp_dict[feat] = []

    # sort the indeces of the features in the vector
    feats_sorted = []
    for feat in feat_order:
        feat_imp = []
        start = feat_order[feat][0]
        end = feat_order[feat][1]
        feat_tuple = (start, end, feat)
        feats_sorted += [feat_tuple]
    feats_sorted.sort()

    # list to hold mean importance value of each index in the vector
    m_index_importance = [0] * vector_length
    sd_index_importance = [0] * vector_length

    # matrix to hold index importance of each index in the vector for all random forest iterations
    f_index_importance = [[]] * vector_length
    for i in range(len(f_index_importance)):
        f_index_importance[i] = []

    ########################################### initiate performance metrics lists
    f_accuracy, f_precision, f_specificity, f_mcc, f_recall, f_log_loss, f_f1_score, f_auc = \
        [], [], [], [], [], [], [], []

    ########################################### initiate TRP and FPR lists for ROC curve
    tpr_list = []
    fpr_list = []
    threshold_list = []
    fpr_mean = []
    tpr_mean = []

    # lists for PR curve
    recall_pr_list = []
    precision_pr_list = []
    average_precision_pr_list = []

    total_test_label_lists = []
    total_prob_label_lists = []

    ########################################### initiate lists for train and test sets during cross-validation
    train_set_iterations = []
    test_set_iterations = []

    ########################################### initiate list for maximum_depths of each tree
    max_tree_depth = []

    ########################################### initiate list for pairs
    f_gene_pairs_list = []
    f_gene_pairs = {}

    ######################################### start iterating the dictionaries
    for trial in range(0, len(set_dicts)):

        # define the current dictionary
        set_dict = set_dicts[trial]

        iterations_out_file.write('RANDOM FOREST' + str(trial + 1) + '\n')

        ################################### creating initialized values for Random Forest predictor

        # prepare datasets for random forest
        instance_set = []  # list with negative and positive instances
        label_set = []  # list with class label
        instance_keys = []  # list for instance IDs as they are added
        p_pair_set = []  # list with tuples of gene pairs for positive instances
        n_pair_set = []  # list with tuples of gene pairs for negative instances


        set_keys = set_dict.keys()
        set_keys = sorted_strings(set_keys)
        test_indeces = []
        test_keys = []

        # append dictionary information in the lists
        for key in set_keys:
            instance_set += [set_dict[key]['comb']]
            instance_keys += [key]
            label_set += [set_dict[key]['class']]
            if 'neg' in key:
                n_pair_set += [set_dict[key]['pair']]
            else:
                p_pair_set += [set_dict[key]['pair']]

        instance_set_array = np.array(instance_set)

        # define pair sets
        final_neg_pairs_list = set()
        for pair in n_pair_set:
            final_neg_pairs_list.add(tuple(pair))
        for pair in p_pair_set:
            final_neg_pairs_list.add(tuple(pair))


        ################### create cross-validation groups (same order with instances)
        instance_groups = stratified_gene_pair_groups(set_dict, instance_keys)

        # create data split
        logo = LeaveOneGroupOut()
        data_split = []
        for train_index, test_index in logo.split(instance_set, label_set, groups=instance_groups):
            data_split += [[train_index, test_index]]

        print('##### Building Random Forest predictor for set %d with the LeaveOnePairOut cross-validation procedure (%d iterations)' % \
        (trial + 1, len(data_split)))


        ################################ create overall sets of iterations for statistics

        sgp_overall_test_set = []
        sgp_overall_test_label = []
        sgp_overall_pred_label = []
        sgp_overall_prob_label = []


        ################################ create feature importance statistics for Random Forest
        cv_feat_imp_dict = {}
        for s, e, f in feats_sorted:
            cv_feat_imp_dict[f] = []

        # matrix to hold importance of each index for all cross validations
        cv_index_importance = [[]] * instance_set_array.shape[1]
        for i in range(len(cv_index_importance)):
            cv_index_importance[i] = []

        # create statistics for the specific Random Forest iteration
        cv_accuracy, cv_precision, cv_recall, cv_specificity, cv_mcc, cv_log_loss, cv_f1_score, cv_auc = \
            [], [], [], [], [], [], [], []

        # lists to hold random forests  and train/test sets for all cross-validations
        rf_models_validation = []
        train_set_validation = []
        test_set_validation = []

        ########################################### train the Random Forest predictor
        for train_index, test_index in data_split:
            test_indeces += [list(test_index)]
            test_keys += [instance_keys[i] for i in test_index]
            train_set = [instance_set[i] for i in train_index]
            test_set = [instance_set[i] for i in test_index]
            train_label = [label_set[i] for i in train_index]
            test_label = [label_set[i] for i in test_index]

            # append sets in the validation sets
            train_set_validation += [list(train_set)]
            test_set_validation += [list(test_set)]

            # initiate the Random Forest instance
            clf = RandomForestClassifier(n_estimators=100, max_depth=10)
            new_clf = clf.fit(train_set, train_label)

            # save rf_model
            rf_models_validation += [new_clf]

            # take maximum tree depth
            max_tree_depth += [max([estimator.tree_.max_depth for estimator in new_clf.estimators_])]

            # see probabilities of predicted classes for test dataset
            prob_classes = new_clf.predict_proba(test_set)  # first element for 0,second for 1

            # see predicted classes of test dataset
            pred_classes = predict_classes(prob_classes, float(parser['class_threshold']))

            # create overall sets for stratified gene pairs cross validation
            test_label = np.array(test_label)
            for el in test_set:
                sgp_overall_test_set += [el]
            for el in pred_classes:
                sgp_overall_pred_label += [el]
            for el in prob_classes:
                sgp_overall_prob_label += [el]
            for el in test_label:
                sgp_overall_test_label += [el]

            ############################### calculate feature importance and accuracy
            feat_importance = new_clf.feature_importances_

            std = np.std([tree.feature_importances_ for tree in new_clf.estimators_], axis=0)

            # append accuracy
            accuracy = new_clf.score(test_set, test_label)
            cv_accuracy += [accuracy]


            # for each feature
            for s, e, f in feats_sorted:
                # take the importance of its indeces
                importance_list = feat_importance[s:e + 1]

                # calculate the sum
                importance_sum = 0
                for el in importance_list:
                    importance_sum += el
                cv_feat_imp_dict[f] += [importance_sum]

            # take the importance of each index in the vector
            for i in range(len(feat_importance)):
                cv_index_importance[i] += [feat_importance[i]]


        ###################################### calculate statistics for cross-validation iterations
        # add train/sets for each set
        train_set_iterations += [train_set_validation]
        test_set_iterations += [test_set_validation]

        # if stratified gene pairs cross validation is selected
        sgp_overall_test_set = np.array(sgp_overall_test_set)
        sgp_overall_pred_label = np.array(sgp_overall_pred_label)
        sgp_overall_prob_label = np.array(sgp_overall_prob_label)


        # append them in the total test and probability labels
        if len(sgp_overall_test_label) == len(sgp_overall_prob_label):

            if len(sgp_overall_test_label) == (combs_num):

                for element in range(0,len(sgp_overall_test_label)):
                    sgp_overall_test_label[element] = np.int32(sgp_overall_test_label[element])

                total_test_label_lists += [list(sgp_overall_test_label)]
                total_prob_label_lists += [list(sgp_overall_prob_label[:, 1])]

        # create confusion matrix
        confusion_matrix = metrics.confusion_matrix(sgp_overall_test_label, sgp_overall_pred_label)
        tp, fp, fn, tn = confusion_matrix[0][0], confusion_matrix[0][1], confusion_matrix[1][0], \
                         confusion_matrix[1][1]


        # calculate statistics
        cv_precision = metrics.precision_score(sgp_overall_test_label, sgp_overall_pred_label)
        cv_recall = metrics.recall_score(sgp_overall_test_label, sgp_overall_pred_label)
        cv_log_loss = metrics.log_loss(sgp_overall_test_label, sgp_overall_prob_label)
        cv_auc = metrics.roc_auc_score(sgp_overall_test_label, sgp_overall_pred_label)
        cv_f1_score = (2 * cv_precision * cv_recall) / (cv_precision + cv_recall)
        cv_specificity = tn / (tn + fp)

        # Matthews correlation coefficient
        cv_mcc = (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

        # turn to lists
        cv_precision = [cv_precision]
        cv_recall = [cv_recall]
        cv_log_loss = [cv_log_loss]
        cv_auc = [cv_auc]
        cv_f1_score = [cv_f1_score]
        cv_specificity = [cv_specificity]
        cv_mcc = [cv_mcc]

        # save TPR and FPR rates for ROC curve
        # find false positive, true positive rates and append in the lists
        fpr_l, tpr_l, thresholds = metrics.roc_curve(sgp_overall_test_label, sgp_overall_prob_label[:, 1],
                                                     pos_label=1)
        fpr_list += [fpr_l]
        tpr_list += [tpr_l]
        threshold_list += [thresholds]

        # save precision and recall rates for PR curve
        precision_l, recall_l, _ = precision_recall_curve(sgp_overall_test_label, sgp_overall_prob_label[:, 1])
        average_precision_l = average_precision_score(sgp_overall_test_label, sgp_overall_prob_label[:, 1])

        recall_pr_list += [recall_l]
        precision_pr_list += [precision_l]
        average_precision_pr_list += [average_precision_l]

        ######################################################################

        # find mean feature significance
        iterations_out_file.write('Feature significance\n')
        cv_feat_importance_list = []
        for feat in cv_feat_imp_dict:
            mean_importance = mean(cv_feat_imp_dict[feat])
            feat_imp_dict[feat] += [mean_importance]
            cv_feat_importance_list += [(mean_importance, feat)]
        cv_feat_importance_list.sort()
        # write to file
        for imp, f in cv_feat_importance_list:
            iterations_out_file.write(f + '\t' + str(imp) + '\n')

        # find mean importance of each index in the vector
        for i in range(len(cv_index_importance)):
            f_index_importance[i] += [mean(cv_index_importance[i])]

        # find mean values of statistics
        f_accuracy += [mean(cv_accuracy)]
        f_precision += [mean(cv_precision)]
        f_recall += [mean(cv_recall)]
        f_specificity += [mean(cv_specificity)]
        f_mcc += [mean(cv_mcc)]
        f_log_loss += [mean(cv_log_loss)]
        f_f1_score += [mean(cv_f1_score)]
        f_auc += [mean(cv_auc)]

        # write statistics to file
        iterations_out_file.write('\nStatistics\n')
        iterations_out_file.write('Accuracy: ' + str(mean(cv_accuracy)) + '\n')

        iterations_out_file.write('Precision: ' + str(mean(cv_precision)) + '\n')
        iterations_out_file.write('Recall: ' + str(mean(cv_recall)) + '\n')
        iterations_out_file.write('Specificity: ' + str(mean(cv_specificity)) + '\n')
        iterations_out_file.write('Matthews Correlation Coefficient: ' + str(mean(cv_mcc)) + '\n')
        iterations_out_file.write('Log loss: ' + str(mean(cv_log_loss)) + '\n')
        iterations_out_file.write('F1 score: ' + str(mean(cv_f1_score)) + '\n')
        iterations_out_file.write('AUC: ' + str(mean(cv_auc)) + '\n')


        # get unique gene pair information for negative dataset for this tree
        cv_gene_pairs = {}  # dictionary only for this tree
        for geneA, geneB in n_pair_set:
            pair_str1 = '-'.join([geneA, geneB])
            pair_str2 = '-'.join([geneB, geneA])

            # append to the dictionary for this tree
            if pair_str1 not in cv_gene_pairs and pair_str2 not in cv_gene_pairs:
                cv_gene_pairs[pair_str1] = 1
            elif pair_str1 in cv_gene_pairs:
                cv_gene_pairs[pair_str1] += 1
            elif pair_str2 in cv_gene_pairs:
                cv_gene_pairs[pair_str2] += 1

            # append to the final dictionary
            if pair_str1 not in f_gene_pairs and pair_str2 not in f_gene_pairs:
                f_gene_pairs[pair_str1] = 1
            elif pair_str1 in f_gene_pairs:
                f_gene_pairs[pair_str1] += 1
            elif pair_str2 in f_gene_pairs:
                f_gene_pairs[pair_str2] += 1

                # append the length of the list in the final list
        f_gene_pairs_list += [len(cv_gene_pairs)]


        # write to the file
        iterations_out_file.write(
            'Number of unique pairs in the negative file: ' + str(len(cv_gene_pairs.keys())) + '\n\n')


    ################################### calculate statistics for all iterations

    summary_out_file.write('####### Random Forest Feature Mean Importance (descending order) \n')
    summary_out_file.write('Feature\tMean\tStandard deviation\n')

    ######################## write fpr, tpr, thresholds
    json_fpr_list = []
    json_tpr_list = []
    json_threshold_list = []

    for el in fpr_list:
        json_fpr_list += [list(el)]
    for el in tpr_list:
        json_tpr_list += [list(el)]
    for el in threshold_list:
        json_threshold_list += [list(el)]

    final_metrics_dict = {'fpr': json_fpr_list, 'tpr': json_tpr_list, 'roc_thresholds': json_threshold_list}

    ######################## write precision,recall
    json_recall_list = []
    json_precision_list = []
    json_precision_average_list = []

    for el in recall_pr_list:
        json_recall_list += [list(el)]
    for el in precision_pr_list:
        json_precision_list += [list(el)]
    for el in average_precision_pr_list:
        json_precision_average_list += [el]

    final_metrics_dict.update({'recall': json_recall_list, 'precision': json_precision_list,
                               'precision_average': json_precision_average_list})

    ######################## write true test labels, probability labels
    total_test_label_lists = list(total_test_label_lists)
    total_prob_label_lists = list(total_prob_label_lists)
    for i in range(len(total_test_label_lists)):
        total_test_label_lists[i] = list(total_test_label_lists[i])
        total_prob_label_lists[i] = list(total_prob_label_lists[i])


    final_metrics_dict.update({'prob_labels': total_prob_label_lists})
    final_metrics_dict.update({'test_sets': test_set_validation, 'train_sets': train_set_validation})


    with open(output_dir + parser['filename'] + '_metrics_dict.json', 'w') as fpr_json:
        json.dump(final_metrics_dict, fpr_json)


    # calculate mean importance of each index in the vector
    for i in range(len(f_index_importance)):
        m_index_importance[i] = mean(f_index_importance[i])
        sd_index_importance[i] = sd(f_index_importance[i], m_index_importance[i])

    m_index_importance = np.array(m_index_importance)
    m_indices = np.argsort(m_index_importance)[::-1]

    # calculate means and sd of statistic measures (in overall)
    mean_f_accuracy = mean(f_accuracy)
    sd_f_accuracy = sd(f_accuracy, mean_f_accuracy)

    mean_f_precision = mean(f_precision)
    sd_f_precision = sd(f_precision, mean_f_precision)
    mean_f_recall = mean(f_recall)
    sd_f_recall = sd(f_recall, mean_f_recall)
    mean_f_specificity = mean(f_specificity)
    sd_f_specificity = sd(f_specificity, mean_f_specificity)
    mean_f_mcc = mean(f_mcc)
    sd_f_mcc = sd(f_mcc, mean_f_mcc)
    mean_f_log_loss = mean(f_log_loss)
    sd_f_log_loss = sd(f_log_loss, mean_f_log_loss)
    mean_f_f1_score = mean(f_f1_score)
    sd_f_f1_score = sd(f_f1_score, mean_f_f1_score)
    mean_f_auc = mean(f_auc)
    sd_f_auc = sd(f_auc, mean_f_auc)

    ################################################# count votes among RFs

    #### collect predictions and probabilities for voting



    # calculate statistics on negative gene pairs
    mean_gene_pair_num = mean(f_gene_pairs_list)
    sd_gene_pair_num = sd(f_gene_pairs_list, mean_gene_pair_num)
    # find number of appearance for each negative gene pair
    f_gene_pairs_count_list = []
    for key in f_gene_pairs:
        f_gene_pairs_count_list += [(f_gene_pairs[key], key)]
    f_gene_pairs_count_list.sort()


    # calculate feature importance,write to file
    feat_importance_list = []
    for feat in feat_imp_dict:
        f_mean_importance = mean(feat_imp_dict[feat])
        f_sd_importance = sd(feat_imp_dict[feat], f_mean_importance)
        feat_importance_list += [(f_mean_importance, f_sd_importance, feat)]

    feat_importance_list.sort(reverse=True)

    # write overall feature importance in the file
    for m, sd, f in feat_importance_list:
        summary_out_file.write(f + '\t' + str(m) + '\t' + str(sd) + '\n')

    # write mean importance of each index in the vector in the file
    summary_out_file.write('\n####### Importance of each vector index, starting from element 1 (descending order) \n')
    summary_out_file.write('Rank\tVector_index\tImportance\tStandard_deviation\tFeature_name\tFeature_index\n')
    for f in range(instance_set_array.shape[1]):
        summary_out_file.write(
            str(f + 1) + '\t' + str(m_indices[f] + 1) + '\t' + str(m_index_importance[m_indices[f]]) + '\t' + str(
                sd_index_importance[m_indices[f]]))
        for s, e, g in feats_sorted:
            if m_indices[f] in range(s, e + 1):
                summary_out_file.write('\t' + g + '\t' + str(range(s, e + 1).index(m_indices[f]) + 1) + '\n')
                break

    ################## plot the feature importance ####################
    # create three different arrays for mean, std and feature names
    bar_mean = []
    bar_std = []
    bar_feats = []

    for m, sd, f in feat_importance_list:
        bar_mean += [m]
        bar_std += [sd]
        bar_feats += [f]

    final_feat_dict = [{'mean': bar_mean, 'std': bar_std, 'feats': bar_feats}, feat_imp_dict]

    with open(output_dir + parser['filename'] + '_feat_dict.json', 'w') as feat_json:
        json.dump(final_feat_dict, feat_json)

    # write statistics to the file

    summary_out_file.write('\n####### Number of features\n')
    summary_out_file.write('Number of distinct features: ' + str(len(feat_imp_dict.keys())) + '\n')
    summary_out_file.write('Number of vector features: ' + str(len(instance_set[0])) + '\n')
    summary_out_file.write('Maximum tree depth: ' + str(max(max_tree_depth)) + '\n')

    summary_out_file.write('\n####### Random Forest Statistics (performance among all 500 RFs) \n')
    summary_out_file.write('Measure\tMean\tStandard deviation\n')
    summary_out_file.write('Accuracy\t' + str(mean_f_accuracy) + '\t' + str(sd_f_accuracy) + '\n')
    summary_out_file.write('Precision\t' + str(mean_f_precision) + '\t' + str(sd_f_precision) + '\n')
    summary_out_file.write('Recall\t' + str(mean_f_recall) + '\t' + str(sd_f_recall) + '\n')
    summary_out_file.write('Specificity\t' + str(mean_f_specificity) + '\t' + str(sd_f_specificity) + '\n')
    summary_out_file.write('Matthews Correlation Coefficient\t' + str(mean_f_mcc) + '\t' + str(sd_f_mcc) + '\n')
    summary_out_file.write('Log loss\t' + str(mean_f_log_loss) + '\t' + str(sd_f_log_loss) + '\n')
    summary_out_file.write('F1 score\t' + str(mean_f_f1_score) + '\t' + str(sd_f_f1_score) + '\n')
    summary_out_file.write('AUC\t' + str(mean_f_auc) + '\t' + str(sd_f_auc) + '\n')

    summary_out_file.write('\n####### Unique Feature statistics on 1KGP instances (among all 500 RFs) \n')
    summary_out_file.write('Feature\tMean\tStandard deviation\n')
    summary_out_file.write('Unique 1KGP gene pairs\t' + str(mean_gene_pair_num) + '\t' + str(sd_gene_pair_num) + '\n')


    summary_out_file.write('\n####### Gene pairs on 1KGP instances\n')
    summary_out_file.write('Gene Pair\tTotal count number\tMean (over all iterations)\n')
    for gene_combi in f_gene_pairs_count_list:
        num = gene_combi[0]
        comb = gene_combi[1]
        summary_out_file.write(comb + '\t' + str(num) + '\t' + str(num / len(set_dicts)) + '\n')

    iterations_out_file.close()































