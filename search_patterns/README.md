#SearchIt1.0 - January 2014 
Author: Francesco Musacchia

This script has the purpose to search elements in specific fields of an Annocript ouput

First thing is to open Annocript output and patterns file given in input
The patterns file is stored in a list of string
The script will read the output and search the patterns for each line. Depending from the patterns found,
to each line is assigned a score. These scores are stored in an array [lines,patterns]. 
The script goes through the matrix to search for K elements with highest scores. K is given in input

USAGE: to use the script please call it giving in input:
     - the filtered output of Annocript 
     - the file with patterns:
           The file with patterns has the first line with the fields where you want to search (separated by ';')
            in second and other lines we have the searches: for each line you should write: words to search = importance
           three levels of importance are given: H: high (1); M: medium (0.5); L: low (0.1). You can use all the words you want
             An example of patterns.txt file:
                   DescriptionSP;DescriptionTR
                   actin=H;cellular=M
                   transposon=L;binding=H;neural=H
                   molecular=L

     - the number of rows to write for each search

OUTPUT: Returns p files (p is the number of searches) with the rows of Annocript where the patterns are found.

EXAMPLE USAGE: python annosearch.py your_filtered_ann_out.txt patterns.txt 100
