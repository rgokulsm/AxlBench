'''
Created on Apr 10, 2011

@author: Hadi Esmaeilzadeh <hadianeh@cs.washington.edu>
'''

import png
import sys
import csv

def rgbq(file1, file2):
    csvReader = csv.reader(open(file1, 'r'), delimiter=',', quotechar='"')   
    csvReader2 = csv.reader(open(file2, 'r'), delimiter=',', quotechar='"')   


    i = 0
    j = 0
    pixels = []
    pixels2= []
    width = 0
    height = 0
    meta = {}
    for row,row2 in csvReader,csvReader2:
        if (i == 0):
            width = int(row[0])
            height = int(row[1])
            print(width, height)
        elif (i == height + 1):
            meta = row[0]
            print(meta)
            break
        else: 
             row = [int(e) for e in row]
	     row2 = [int(f) for f in row2]
             for e,f in row,row2:
		print (e,f)
        pass
        
        i = i + 1
    pass


#    for row in csvReader2:
#        if (i == 0):
#            width = int(row[0])
#            height = int(row[1])
#            print(width, height)
#        elif (i == height + 1):
#            meta = row[0]
#            print(meta)
#            break
#        else: 
##            row = [int(e) for e in row]
#             for e in row:
#		pixels2.append(int(e))
#        pass
#        
#        i = i + 1
#    pass
#
#
#for f, b in zip(pixels, pixels2):
#    print(f, b)
#
#
    
#    print(width, height, meta)
#    return(width, height, pixels, meta)
pass



if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print('Error: Oops! Too few arguments!')
        print('Usage: ' + sys.argv[0] + ' OPERATION INPUT_FILE OUTPUT_FILE')
        exit(-1)
    pass

    input = str(sys.argv[1])
    precise = str(sys.argv[2])

    rgbq(input, precise)

    
    exit(0)

pass
