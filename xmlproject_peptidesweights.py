'''
Rebecca Fuchs
Project 2-MS/MS Viewer!!!!
This code will take command line input

command line: python xmlproject_peptidesweights.py mzxml.file scan.number peptide.seq
python xmlproject_peptidesweights.py 17mix_test2.mzxml 1298 TYDSYLGDDYVR

It will print the matching peaks on the command line and will output a graph showing
the peaks that match the b/y ions caluculated by the code
'''

#command line code
import sys
try:
    input1=sys.argv[1]
    input2=sys.argv[2]
    input3=sys.argv[3]
except IndexError:
    print ("You forgot to provide mzxml.file OR scan.number OR peptide.seq")
    sys.exit(1)

#######################################################getting b y ion values

# a dictionary of the  molecular weights of amino acids 
aaweights = {
	'A' : 71.04,  # alanine
	'R' : 156.10, # arginine
	'D' : 115.03, # aspartic acid
	'N' : 114.04, # asparagine
	'C' : 103.15, # cysteine
	'E' : 129.04, # glutamic acid
	'Q' : 128.06, # glutamine
	'G' : 57.05,  # glycine
	'H' : 137.14, # histidine
	'I' : 113.16, # isoleucine
	'L' : 113.16, # leucine
	'K' : 128.10, # lysine
	'M' : 131.04, # methionine
	'F' : 147.1, # phenylalanine
	'P' : 97.12,  # proline
	'S' : 87.08,  # serine
	'T' : 101.11, # threonine
	'W' : 186.1, # tryptophan
	'Y' : 163.06, # tyrosine
	'V' : 99.07   # valine
}


#below is peptide input, must be made a LIST!!
peptide_try=list(input3)
print("The peptide sequence you input", peptide_try)
bion_dict={}
yion_dict={}

plength=len(peptide_try)
#print(plength)

#code to get b/y ion weights, loop through each segment and get mass for each aa
for i in range(1,plength+1):
    bweight=0
    yweight=0
    bion=peptide_try[0:i]
    yion=peptide_try[i-1:plength+1]
    #print (i)
    for aab in bion:
        bweight=aaweights[aab]+bweight
    bion_dict['%db'%i]=bweight+1

    for aay in yion:
        yweight=aaweights[aay]+yweight
    yion_dict['%dy'%i]=yweight+19
    
        



#make python dictionaries into lists
bion_list=list(bion_dict.values())
yion_list=list(yion_dict.values())



bion_len=len(bion_list)
bionints=[500]*bion_len #intensity of b/y ions needs to be large enough to see on graph
yion_len=len(yion_list)
yionints=[500]*yion_len



############################################### dealing with the xml file
import xml.etree.ElementTree as ET

try:
    tree = ET.parse(input1)
    root = tree.getroot()
except FileNotFoundError:
    print ("File not found, misspelled or in wrong directory")
    sys.exit(1)
#print(root.tag)
l=[]


#xml parsing
d={}

scannum=input2 #scan number
stringnum=str(scannum)

for child in root:
    if 'num' in child.attrib:
        if child.attrib['num'] ==str(scannum):
            peaks_binary=child.attrib
            #print (child.attrib)
            for peaks in child:
                #print (peaks.text)
                peaks_text_binary=peaks.text
                       
############################################decoding

from base64 import b64decode
from array import array

   # peaks elt is the XML element corresponding to the peaks list
try:
    peaks = array('f',b64decode(peaks_text_binary))

    if sys.byteorder != 'big':
        peaks.byteswap()
except NameError:
    print("Scan number does not exist:"+stringnum)
    sys.exit(1)
#print("mzs are below---")
mzs = peaks[::2]
#print(mzs)

#print('ints are below---')
ints = peaks[1::2]
#print(ints)

####################################compare mzs to bion dict values
bion_dict_matches={}
mzs_list=[]
#make mzs into list
for i in range(1,(len(ints)+1)):
    mzs_list.append(peaks[2*(i-1)])
#print (mzs_list)   

#if mzs is within 0.2 of the b ion it will be added to the match dictionary
for key, value in bion_dict.items():
    #print(value)
    for mz in mzs_list:
        #print(mz)
        if (mz-value) > -0.2 and (mz-value)<0.2:
          
            #print (mz,value)
            bion_dict_matches[key]=value
         
        
###################################compare mzs to yion dict values
yion_dict_matches={}  
#if mzs is within 0.2 of the y ion it will be added to the match dictionary
for key, value in yion_dict.items():
    #print(value)
    for mz in mzs_list:
        #print(mz)
        if (mz-value) > -0.2 and (mz-value)<0.2:
          
            #print (mz,value)
            yion_dict_matches[key]=value
print ("Here are the b & y ion matches, they are labeled in the figure too")
print (bion_dict_matches)
print (yion_dict_matches)

'''
some notes  :)
#here are the mzs values, ints and mzs are the same length so im guna put in dict?
for i in range(1,(len(ints)+1)):
    print(peaks[2*(i-1)])

print("=======")
#here are the intensitys values
for i in range(1,(len(ints)+1)):
    print(peaks[2*(i-1)+1])

?peaks[2*(i-1)]   is the m/z value of the ith peak (i=1,2,3,4,...)
?peaks[2*(i-1)+1] is the intensity value of the ith peak (i=1,2,3,4,...)
!!!mzs[i-1] is the m/z value of the ith peak
!!!ints[i-1] is the intensity value of the ith peak
'''


#######################################plotting time 


from pylab import *
import matplotlib.pyplot as plt
import numpy as np

#peaks from scan
(markers, stemlines, baseline) = plt.stem(mzs, ints, markerfmt=" ")
plt.setp(stemlines, linestyle="-", color="blue", linewidth=1)
plt.setp(baseline, visible=False)

#peak from bion
(markers, stemlines, baseline) = plt.stem(bion_list, bionints, markerfmt=" ")
plt.setp(stemlines, linestyle="-", color="red", linewidth=0.5)
plt.setp(baseline, visible=False)

#label for matches
for key, value in bion_dict_matches.items():
    plt.text(value, 501, key,fontsize=5)

#peak from yion   
(markers, stemlines, baseline) = plt.stem(yion_list,yionints, markerfmt=" ")
plt.setp(stemlines, linestyle="-", color="red", linewidth=0.5)
plt.setp(baseline, visible=False)

#label for matches
for key, value in yion_dict_matches.items():
    plt.text(value, 501, key,fontsize=5)

#put the peptide on the figure
#plt.text(right,bottom, peptide_try, fontsize=5, ha="right", color='.5', horizontalalignment='center',verticalalignment='top', transform=ax.transAxes)
#doesnt look right for all figures :(
    
#labels 
plt.xlabel("m/z")
plt.ylabel("relative abundance")
plt.title('Peptide mass spec for scan:'+stringnum )
plt.legend (( 'mass spec shown in blue', 'b/y bions show in red'), loc='upper center', shadow=True)
plt.show()

