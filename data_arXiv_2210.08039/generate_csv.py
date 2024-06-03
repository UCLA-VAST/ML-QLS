import os
import csv

myPath = "/home/wanhsuan/mOLSQ/sabre_log/"

files = [f for f in os.listdir(myPath)]
# print(files)


with open('sabre_results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    field = ["file name", "gate num", "swap", "depth", "time"]
    
    writer.writerow(field)

    for file in files:
        # print(file)
        f = open(myPath+file, "r")
        lines = f.readlines()
        row = [file, None, None, None, None]
        for line in lines:
            line = line.split()
            if line[0] == "gate":
                row[1] = int(line[2])
            elif line[0] == "Time:":
                row[4] = float(line[1])
            elif line[0] == "count_swap:":
                row[2] = int(line[1])
            elif line[0] == "depth:":
                row[3] = int(line[1])
        writer.writerow(row)

# gate num: 1342
# Time:  0.10282554803416133
# count_swap:  165
# depth: 130