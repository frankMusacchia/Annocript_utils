#AnnoSearch1.0 - January 2014 
#Author: Francesco Musacchia
#This script has the purpose to search elements in specific fields of an Annocript ouput

# First thing is to open Annocript output and patterns file given in input
# The patterns file is stored in a list of string
# The script will read the output and search the patterns for each line. Depending from the patterns found,
# to each line is assigned a score. These scores are stored in an array [lines,patterns]. 
# The script goes through the matrix to search for K elements with highest scores. K is given in input

#USAGE: to use the script please call it giving in input:
#       - the filtered output of Annocript 
#       - the file with patterns:
#             The file with patterns has the first line with the fields where you want to search (separated by ';')
#              in second and other lines we have the searches: for each line you should write: words to search = importance
#             three levels of importance are given: H: high (1); M: medium (0.5); L: low (0.1). You can use all the words you want
#             An example of patterns.txt file:
#                   DescriptionSP;DescriptionTR
#                   actin=H;cellular=M
#                   transposon=L;binding=H;neural=H
#                   molecular=L
#
#       - the number of rows to write for each search
#
#OUTPUT: Returns p files (p is the number of searches) with the rows of Annocript where the patterns are found.


#IMPORT
import sys, getopt
import numpy
import linecache #To access the annocript out file at a chose line
import re #Import regexp library to search whole words
import os.path #To access files



#GLOBAL
ann_out_path = ""#Annocript output file
patterns_file_path = "" #Patterns file
scores = [] #Matrix of scores for each pattern and row
tot_patterns = 0 #Total of patterns
ann_out_lines = [] #the rows
ann_out_first_line = ""#We need the first line of the output file to give headers
patterns_names = [] #List of pattern names to give names at the files
k = 0

#Here I check the input parameters to the program
def check_args():
  global ann_out_path
  global patterns_file_path
  global k

  if int(len(sys.argv)) < 3:
    print "Usage: %s annocript_out_file patterns_file results_number" %(sys.argv[0])
    sys.exit(0)
    
  if sys.argv[1] and os.path.isfile(sys.argv[1]):
    ann_out_path = sys.argv[1] #Annocript output file
  else:
    print "You should give the path to an Annocript output! Exiting..."
    sys.exit(0)
    
  if sys.argv[2] and os.path.isfile(sys.argv[2]):
    patterns_file_path = sys.argv[2] #Patterns file
  else:
    print "SearchIt1.0 wants a file with pattern and fields to search. Please create it!"
    sys.exit(0)
  
  #Checking if k is given
  if int(len(sys.argv)) == 4:
    k = int(sys.argv[3])
  else:
    k = 10 
    print "Number of results wanted has not been given. Default is %d" %(k)
   


#STORAGE OF THE PATTERNS and THE FIELDS TO SEARCH
#Here we create three lists: 1. the list of the rows from the pattern file; 2. the list of fields to search
# 3. a list with names taken from patterns to give names to files
def store_patterns (patterns_file_path):
  fields = []
  patterns_list = []
  patterns_file = open (patterns_file_path,"r")
  global patterns_names
  
  
  #Extract a list with the field to search
  fields_line = patterns_file.readline()
  fields_line = fields_line.rstrip('\n')
  fields = fields_line.split(";")
  
  #extract a list with the rows of the pattern file
  lines = patterns_file.readlines()
  for line in lines:
    line = line.rstrip()
    patterns_list.append(line)
    
    #Extract the first name to give a name to the search out file
    patterns_names.append( (line.split(";")[0]).split("=")[0] )
          
  return fields,patterns_list


#DEFINITION OF SCORES
#This permits to work like with a SWITCH.. CASE construct 
def get_score (letter):
  score = 0
  if letter == "H":
    score = 1
  elif letter == "M":
    score = 0.5
  elif letter == "L":
    score = 0.1
  return score
  
#Find whole word using a regex
#The re.search function search string1 in string2 by using \b to find whole word and ignoring the case 
def whole_string_found(string1, string2):
  #return re.search(r'\b%s\b' % (re.escape(string1)), string2,re.IGNORECASE) is not None
  return re.search(r'\b%s\b' % (string1), string2,re.IGNORECASE) is not None

####################
#This function creates an array with scores. For each line and each pattern there'll be a score.
#In this way we should only sort the cols separately and fetch the first k indexes to get the wanted result    
def assign_scores ():
  #Redefinition of global variables: in Python you should always redefine locally and start with 'global'
  global tot_patterns
  global scores
  global ann_out_lines
  global ann_out_first_line
  
  logFile = open("log.txt","w")
  
  #Let's open the Annocript output
  ann_out_file = open(ann_out_path,"r")
  
  #Read any line
  ann_out_first_line = ann_out_file.readline()
  ann_out_lines = ann_out_file.readlines()
  
  #Assign index to header fields
  field_indexes = {}
  field_counter = 0
  ann_out_fields = ann_out_first_line.split("\t")
  for ann_out_field in ann_out_fields:
    field_indexes[ann_out_field] = field_counter
    field_counter+=1
  
  #Store the patterns and the fields to search in a couple of lists
  fields_to_search,patterns = store_patterns(patterns_file_path)
  
  print '[%s]' % ', '.join(map(str, fields_to_search))
  
  #Initialize counters
  ann_out_line_num = 1#Starts from 1 without the header
  pattern_num = 0
  print "Annocript output lines: %d and patterns to search: %d " %(len(ann_out_lines),len(patterns))
  
  #Initialize an array that for each pattern and each row assigns a score
  scores = [[0 for i in range(len(ann_out_lines))] for j in range(len(patterns))]
  
  tot_patterns = len(patterns)
  
  #Loop for each line of the ann out
  for ann_out_line in ann_out_lines:
    ann_out_line_fields = ann_out_line.split("\t")
    
    #Create a list with the fields to be searched
    fields = []
    for field_to_search in fields_to_search:
      fields.append(ann_out_line_fields[ field_indexes[field_to_search] ])
    
    #Loop on patterns
    for pattern in patterns:
      #print str(pattern_num)+" pattern: "+pattern+"\n"
      values = pattern.split(";")#Extract patterns from the line
      #Loop on the fields to be searched
      for field in fields:
        #for each pattern
        for value in values:
          val_and_qual = value.split("=")#separate value and score
          val_and_qual[1] = val_and_qual[1].rstrip("\n")
          
          #search the  pattern in the field
          #print "Searching %s in line %d. Quality: %s" %(val_and_qual[0],ann_out_line_num, val_and_qual[1])
          if whole_string_found(val_and_qual[0],field):
            print "Found %s in row %d: %s" %(val_and_qual[0],ann_out_line_num,field)
            logFile.write("Found %s in row %d: %s \n" %(val_and_qual[0],ann_out_line_num,field))
            scores[pattern_num][ann_out_line_num] +=  get_score(val_and_qual[1])
      pattern_num+=1
    ann_out_line_num+=1
    pattern_num = 0
 # print scores
    logFile.close


#Prints on k files (with name "filename_k") the k best rows from the output of Annocript
#given by the scores stored in the 'scores' matrix
def print_best_k (filename,k):
  
  #Initialization of variables
  best = -1
  bestIndex = -1
  best_scores = []
  
  #basename
  out_search_file = filename
  
  #Set the ranges to use for loops
  to_take = range(k)
  rows = range(len(ann_out_lines))
  num_patt = 0
  for pattern in patterns_names:#For each pattern it writes a file
    out_search_file += pattern 
    file = open(out_search_file,"w")
    file.write(ann_out_first_line)
    #Loop for k elements to take
    for i in to_take:
      #For each row we will see the best elements
      for row in rows:
        if scores[num_patt][row] > best:
          bestIndex = row
          best = scores[num_patt][row]
      if best > 0:
        file.write(linecache.getline(ann_out_path, bestIndex+1))
      scores[num_patt][bestIndex] = -100  
      best = -10
    file.close()
    out_search_file = filename
    num_patt+=1
  
    
#MAIN

#Checks if the arguments are given
check_args()

#Assign scores
print ann_out_path+" and "+patterns_file_path
assign_scores()

#Print the best scores for each pattern in some files
print patterns_names
print_best_k("search_",k)
