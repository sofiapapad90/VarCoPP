#!/usr/bin/env python
 
'''

   Author: Sofia Papadimitriou
   
   VarCoPP predicts the pathogenicity of bi-locus variant combinations and provides a
   Classification and Support score for each prediction.

   ** NOTE: This code uses already annotated bi-locus combinations.

   == INPUT == : a .txt TAB-delimited file with the annotated bi-locus combinations, one per line. The first column
   should be the combination ID followed by 11 columns for the corresponding features in the order:
   Flex1, Hydr1, CADD1, CADD2, CADD3, CADD4, HI_A, RecA, HI_B, RecB, Biol_Dist.


   == PROCESS ==: The script uses the already created VarCoPP model (./models/VarCoPP_model.p.gz) to predict the pathogenicity of each bi-locus
   combination in the file.


   == OUTPUT ==: the predictions.txt file with a ranked list of bi-locus combinations based on their
    Support and Classification Scores


   
   Usage examples
   ===========================================================================   
   
   To run VarCoPP.py with the default parameters for the "Dataset_S1.txt" file, you can simply call the
   script:
   $ python VarCoPP.py -v ./validation_data/Dataset_S1.txt

   Make sure that if the file is not in the same directory as VarCoPP.py, you always provide each path.

   You can specify a prefix for the output file like this:
   $ python VarCoPP.py -v ./validation_data/Dataset_S1.txt -f test

   Then, the output file will be named: test_predictions.txt
   

   '''

from __future__ import division
import os.path
import argparse
import tarfile
import pickle
import numpy as np
import gzip


def input_parser():
 
    """ Function to construct and parse input parameters """
 
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-v','--combs_file', 
                        help = 'File with bi-locus combinations')
    parser.add_argument('-m','--model', 
                        help = 'File with trained VarCoPP model (DEFAULT: ./models/VarCoPP_model.p.gz)',
                        default='./models/VarCoPP_model.p.gz')
    parser.add_argument('-f','--filename', 
                        help = 'Specify the prefix name of the output file',
                        default='')

    return vars(parser.parse_args())


def open_file(filename):
    
    ''' 
    
    Function that takes the full name of a file and opens it for reading or 
    prints an error otherwise.
    
    '''

    #try to open input file  
    try:
        if filename.endswith('.tar') or filename.endswith('.gz'):
            files=tarfile.open(filename, 'r') 
            for member in files:       
                input_file=files.extractfile(member)
        else:
            input_file=open(filename)
        return input_file
    
    #if the file does not exist, return an Error
    except IOError:
        print('##### File or directory path %s does not exist'%(filename))
        print('##### VarCoPP is terminated')
        quit()


def predict_classes(class_probs,threshold):
    
    '''Function that takes the probabilities of class assignments for each
    element and the threshold for class 1 and returns an array with class
    labels for each element in the list'''
    
    #list to hold labels
    predicted_classes=[]
    
    #iterate over the elements
    for el in class_probs:
        #check threshold for class 1
        if el[1]>threshold:
            predicted_classes+=[1]
        else:
            predicted_classes+=[0]
    
    predicted_classes=np.array(predicted_classes)
    
    return predicted_classes
    

def majority_vote(class_0,class_1):
    
    ''' 
    Function that takes an array with the true and predicted labels of the instance, 
    among all iterations and provides the majority vote for this instance.

    '''    
    
    #check which class has more trees
    if class_0>=class_1:
        majority_vote=0
    elif class_0<class_1:
        majority_vote=1

    return majority_vote     

def percentage_vote(predicted_probabilities,predicted_classes,true_keys):
    
    ''' Function that takes an array with the true labels of test data and a matrix
    with iterations of predicted classes of the test data (in the same order),
    as well as the assigned probabilities and returns a dictionary with 
    information on how each instance is voted for all iterations
    and a list with the majority votes in the same order as the labels'''
    
    #initiate vote dictionary
    vote_dict={'positive':{},'negative':{}}    
    
    #initiate majority votes
    majority_votes=[]
    majority_median_class1_probs=[]    
    
    #iterate over the elements
    for i in range(len(true_keys)):
        
        all_class_iterations=[] #list to hold all class predictions for this element        
        probability_0=[]
        probability_1=[]    
        true_key=true_keys[i]
       
        #iterate over the iterations
        for j in range(len(predicted_classes)):
            
            #define predicted class
            class_iteration=predicted_classes[j][i]
            all_class_iterations+=[class_iteration] #append all class predictions
            
            #define probabilities
            prob_iteration=predicted_probabilities[j][i]            
                        
            #append probabilities  
            probability_0+=[prob_iteration[0]]
            probability_1+=[prob_iteration[1]]
        
        #count votes        
        num_0=all_class_iterations.count(0)
        num_1=all_class_iterations.count(1)
        
 
        perc_0=(num_0*100)/len(predicted_classes)
        perc_1=(num_1*100)/len(predicted_classes)
        
        #count median probabilities
        probability_0_median=np.median(probability_0)        
        probability_1_median=np.median(probability_1)        
        
        #count maximum probabilities
        max_probability_0=max(probability_0)
        max_probability_1=max(probability_1)        
       
        #find the majority vote for this instance
        majority_vote_class=majority_vote(num_0,num_1) 
    
        majority_votes+=[majority_vote_class]
        majority_median_class1_probs+=[probability_1_median]
        
        #append information       
        if majority_vote_class==1:
            vote_dict['positive'].update({true_key:{'predicted_label':majority_vote_class,
            'count_1':num_1,'count_0':num_0,'perc_0':perc_0,'perc_1':perc_1,'prob_0':probability_0_median,
            'prob_1':probability_1_median,'max_prob_0':max_probability_0,'max_prob_1':max_probability_1}})
            
        else:
            vote_dict['negative'].update({true_key:{'predicted_label':majority_vote_class,
            'count_1':num_1,'count_0':num_0,'perc_0':perc_0,'perc_1':perc_1,'prob_0':probability_0_median,
            'prob_1':probability_1_median,'max_prob_0':max_probability_0,'max_prob_1':max_probability_1}})

    return vote_dict, majority_votes, majority_median_class1_probs


def rank_data(votes_dictionary):
    
    ''' Function that takes a dictionary of information about votes and ranks
    the instances from more-probably disease causing to less probably,
    based first on the median probability and second, on the 
    number of trees agreeing to the classification.
    '''

    #take positive and negative keys
    pos_keys=list(votes_dictionary['positive'].keys())
    neg_keys=list(votes_dictionary['negative'].keys())
    
    #take info on Classification score and Support score
    prob_pos=[]
    prob_neg=[]    
     
    for key in pos_keys:    
            
        prob_pos+=[(votes_dictionary['positive'][key]['perc_1'],                    
                    votes_dictionary['positive'][key]['prob_1'],
                    pos_keys.index(key))]
        
    for key in neg_keys:        
        
        prob_neg+=[(votes_dictionary['negative'][key]['perc_1'],                    
                    votes_dictionary['negative'][key]['prob_1'],                    
                    neg_keys.index(key))]
    
   
    #sort the instances first Support score and then on Classification Score    
    sorted_keys=[]
    
    prob_pos.sort(reverse=True)
    prob_neg.sort(reverse=True)
    
    #create sorted keys
    for el in prob_pos:
        sorted_keys+=[pos_keys[el[2]]]
        
    for el in prob_neg:
        sorted_keys+=[neg_keys[el[2]]]
        
    return sorted_keys,prob_pos,prob_neg,pos_keys,neg_keys
    
def feat_to_vector(al,feat_list):
    
    '''    
    Function that takes a list of features and transforms them into computer
    readable values'    
    '''
 
    #initiate the vector
    vector=[]        
    h=0 #vector index 
   
    #initiate dictionary with the vector start position of each feature
    feat_pos_dict={}
 
    #iterate over the features to vectorize them   
    for i in range(len(feat_list)):    
            
        #vectorize for the aminoacid property change
        if feat_list[i] in ['Flex1','Hydr1']:
            
            #Categorize value or wt
            try:
                if al[i] in ['wt','wild_type']:
                    vector+=[0.0]  
                elif al[i]=='NA':
                    vector+=[vector[-1]]
                else:               
                    vector+=[float(al[i])]   
            except ValueError:
                print('##### Some values of Hydr1 or Flex1 are invalid.'+\
                ' Please use the calc_aa_diff.py to annotate variants'+\
                ' with Hydr and Flex amino acid differences.')
                print('##### Script is terminated')
                quit()
            
            #save vector position of the feature
            feat_name=feat_list[i]
            feat_pos_dict[feat_name]=[h,len(vector)-1]
            h=len(vector)

        #vectorize for CADD
        elif feat_list[i] in ['CADD1','CADD2','CADD3','CADD4']:
            
            #Categorize value or wt  
            try:
                if al[i] in ['wt','wild_type']:
                    vector+=[-3.0]    
                elif al[i]=='NA':
                    vector+=[vector[-1]]
                else:
                    vector+=[float(al[i])]   
            except ValueError:
                print('##### Some CADD raw score values are invalid.'+\
                    ' Please provide a score for CADD for all variants.')
                print('##### Script is terminated')
                quit()
            
            #save vector position of the feature
            feat_name=feat_list[i]
            feat_pos_dict[feat_name]=[h,len(vector)-1]
            h=len(vector)

    ############################ collect information on genes and pair info
    
    #iterate over the gene and pair features
    for i in range(len(feat_list)):  
  
        ############# Haploinsufficiency probability
        if feat_list[i] in ['HI_A','HI_B']:
            
            try:
                if al[i] in ['N/A','NA','NaN','nan']:
                    vector+=[0.19898]
                else:
                    vector+=[float(al[i])]
            except ValueError:
                print('##### Some gene haploinsufficiency values are invalid.'+\
                        ' Please provide a haploinsufficiency value or write'+\
                        ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()
            
            #save vector position of the feature
            feat_pos_dict[feat_list[i]]=[h,len(vector)-1]
            #initiate vector index
            h=len(vector)
        
        ############# Recessiveness probability
        elif feat_list[i] in ['RecA','RecB']:
            
            try:
                if al[i] in ['N/A','NA','NaN','nan']:
                    vector+=[0.12788]
                else:
                    vector+=[float(al[i])]
            except ValueError:
                print('##### Some gene recessiveness values are invalid.'+\
                        ' Please provide a recessiveness value or write'+\
                        ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()
            
            #save vector position of the feature
            feat_pos_dict[feat_list[i]]=[h,len(vector)-1]
            #initiate vector index
            h=len(vector)

        ############ Biological distance
        elif feat_list[i] in ['Biol_Dist']:
            
            try:
                if al[i] in ['N/A','NA','NaN','nan']:
                    vector+=[11.44]
                else:
                    vector+=[float(al[i])]
            except ValueError:
                print('##### Some Biological distance values are invalid.'+\
                        ' Please provide a value or write'+\
                        ' "NA","NaN", or "nan".')
                print('##### Script is terminated')
                quit()
            
            #save vector position of the feature
            feat_pos_dict[feat_list[i]]=[h,len(vector)-1]
            #initiate vector index
            h=len(vector)

    return vector, feat_pos_dict   
   
if __name__=='__main__':         
    
    #Calling the parser function to start parsing the argument
    parser = input_parser()  
    
    ################################################## input checks
    if not os.path.exists(parser['combs_file']):    
        print('##### Input file %s or directory path does not exist'%\
        (parser['combs_file']))
        print('##### VarCoPP is terminated')
        quit()

    ################################################## store combinations
    print('##### Collecting bi-locus combinations')
    
    #initiate lists and dictionaries to hold combinations
    test_set=[]
    test_keys=[]
    test_dict={}
    
    #features list
    feat_list=['Flex1','Hydr1','CADD1','CADD2','CADD3','CADD4','HI_A','RecA']
    feat_list+=['HI_B','RecB','Biol_Dist']
    
    combs_file=open_file(parser['combs_file'])
    combs_header=combs_file.readline().rstrip('\n').rstrip('\r')
    
    #check file
    if len(combs_header.split('\t'))<12:
        print('##### A column is missing from the combinations file header or'+\
        ' the file is not TAB delimited')
        print('##### VarCoPP is terminated')
        combs_file.close()
        quit()
    
    #parse file and collect combinations
    for line in combs_file:
        line_list=line.rstrip('\n').rstrip('\r').split('\t')
        
        if len(line_list)>11:
            
            #create vector
            initial_vector=line_list[1:]
            vector, feat_order =feat_to_vector(initial_vector,feat_list)  
                     
            test_keys+=[line_list[0]]
            test_set+=[vector]        
            
            for i in range(len(vector)):
                vector[i]=str(vector[i])                
            test_dict[line_list[0]]=vector   

    combs_file.close()
   
    ############################################### load Random Forests
    print('##### Getting the model (it takes around 1 1/2 minute)')
    
    try:
        if parser['model'].endswith('.gz'):
            with gzip.open(parser['model'], 'rb') as f:
                model = pickle.load(f)
        else:
           
            with open(parser['model'], 'rb') as f:
                model = pickle.load(f)
    except IOError:
        print('##### VarCoPP trained model could not upload. Make sure that'+\
        ' the model file is in the same directory with VarCoPP.py')
        print('##### VarCoPP is terminated')
        quit()
    
    #initiation of lists for majority vote
    total_test_label_lists=[]   
    total_allclass_prob_lists=[]
    
    ######################################### Predict digenic combinations
    print('##### Predicting... ')
    
    for trial in range(0,len(model)):  
        
        #initiate the Random Forest instance
        trained_model=model[trial][0]
 
        #see probabilities of predicted classes for test dataset 
        #first element for 0, second for 1       
        prob_classes=trained_model.predict_proba(test_set) 
        
        #see predicted classes of test dataset
        pred_classes=predict_classes(prob_classes,0.4891)
 
        #append predicted values in the total test and probability labels
        total_test_label_lists+=[pred_classes]  
        total_allclass_prob_lists+=[prob_classes]
 
    ######################################### calculate majority vote 
     
    vote_dictionary,majority_votes, majority_prob_votes=\
    percentage_vote(total_allclass_prob_lists,total_test_label_lists,
                    test_keys) 
 
    ######################################### rank based on SS and CS
    print('##### Ranking and saving predictions')
    
    ranked_keys,pos_ranked,neg_ranked,pos_keys,neg_keys=\
    rank_data(vote_dictionary)    
    
    #write information in the file
    if parser['filename']=='':
        ranked_file=open('predictions.txt','w')
    else:
        ranked_file=open(parser['filename']+'_predictions.txt','w')
        
    ranked_file.write(combs_header+'\t')
    ranked_file.write('Classification_score\tSupport_score\tPredicted_class\n')
    
    for el in pos_ranked:
        real_key=pos_keys[el[2]]        
        ranked_file.write(real_key+'\t')
        combination='\t'.join(test_dict[real_key])
        ranked_file.write(combination)
        ranked_file.write('\t'+str(el[1])+'\t'+str(el[0])+'\tDisease_causing\n')
        
    for el in neg_ranked:
        real_key=neg_keys[el[2]]  
        ranked_file.write(real_key+'\t')
        combination='\t'.join(test_dict[real_key])
        ranked_file.write(combination)
        ranked_file.write('\t'+str(el[1])+'\t'+str(el[0])+'\tNeutral\n')
    
    ranked_file.close() 