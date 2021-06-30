from csv import DictReader

'''Bismuth Scan Notes COMPLETE - Useful scans
scanindex	date	scan	thickness	energy	beamwaist	fluence	exposure
1	2021-02-17  	11  	23      	1.5 	412 	    0.78	10
2	2021-02-17  	8   	23      	2.5 	412     	1.30	10
3	2021-02-16  	6   	23      	5   	412      	2.60	10
4	2021-02-16  	3   	23      	10  	412     	5.20	10
5	2021-02-17  	13  	23      	15  	412     	7.80	10
6	2021-02-17  	15  	23      	20  	412     	10.4	10
7	2021-02-17  	16  	23      	30  	412     	15.6	10
8	2021-02-21  	5   	14      	10  	412     	5.20	5
9	2021-02-21  	14  	14      	20  	412     	10.4	5
10	2021-02-21  	13  	14        	30  	412     	15.6	5
11	2021-02-21  	12  	14   	    40  	412     	20.8	5
12	2021-02-21  	10  	14      	50  	412     	26.0	5
'''

# filename = 'D:\\Bismuth Project\\New Bismuth Data\\Bismuth Scan Notes COMPLETE - Useful scans.csv'
#
# with open(filename, 'r') as read_obj:
#     dict_reader = DictReader(read_obj)
#     scanlist = list(dict_reader)

filename2 = 'D:\\Bismuth Project\\New Bismuth Data\\Bismuth Scan Notes COMPLETE - Scans by fluence.csv'

with open(filename2, 'r') as read_obj:
    dict_reader = DictReader(read_obj)
    scanlist = list(dict_reader)