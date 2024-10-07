import sys
import os
import re
import pickle

input_file = sys.argv[1]
id = sys.argv[2]
output_dir = os.path.dirname(os.path.abspath(input_file))
hash = {}
global plusminus
plusminus = 2

f = open(input_file, "r")
inputfile = f.read()
for line_num, line in enumerate(inputfile.splitlines(), start=1):
    line = line.strip()
    fields = line.split('\t')
    if len(fields) != 6:
        print(f"Error: Line {line_num} does not have 6 columns")
        print(f"Line content: {line}")
        continue  # Skip lines that don't have 6 fields
    read, ref, read_start, read_end, ref_start, ref_end = fields
    rest = f"{read_start},{read_end},{ref_start},{ref_end}"
    if read not in hash:
        hash[read] = {}
    if ref not in hash[read]:
        hash[read][ref] = {}
    hash[read][ref][int(ref_start)] = rest

readHashc = {}

def junction_looper(r1, r2, rln, junctions, key, key2, cp_key2):
    minr, maxr = (r1, r2) if int(r1) < int(r2) else (r2, r1)
    for each_junction in junctions:
        pl2 = int(each_junction) + plusminus
        mn2 = int(each_junction) - plusminus
        if int(maxr) > pl2 and int(minr) < mn2:
            if 'x' in key2:
                col23 = key2
            else:
                col23 = f"{cp_key2}#{each_junction}"
            
            # Initialize the nested dictionaries if they don't exist
            readHashc.setdefault(key, {}).setdefault(rln, set())
            
            # Add the reference to the set for this read_id and read length
            readHashc[key][rln].add(col23)
            
           # break  

# process hash table
for key in sorted(hash.keys()):
    for key2 in sorted(hash[key].keys()):
        if "x" in key2:
            cp_key2 = key2
        else:
            cp_key2 = re.sub("[,#].+$", "", key2)
        part1, part2 = key2.split("#")
        junctions = part2.split(",")
        entries = []
        for key3 in sorted(hash[key][key2].keys()):
            entries.append(hash[key][key2][key3])
        count = len(entries)
        if count < 2:
            q1, q2, r1, r2 = entries[0].split(",")
            rln = abs(int(r1) - int(r2))
            junction_looper(r1, r2, rln, junctions, key, key2, cp_key2)
            continue
        for x in range(0, count):
            each_entry = entries[x]
            q1, q2, r1, r2 = each_entry.split(",")
            min_each, max_each = (q1, q2) if int(q1) < int(q2) else (q2, q1)
            next_entry = entries[x + 1] if x + 1 < count else ""
            q11, q22, r11, r22 = next_entry.split(",") if next_entry != "" else (None, None, None, None)
            min_next, max_next = (q11, q22) if q11 and int(q11) < int(q22) else (q22, q11) if q22 and int(q22) < int(q11) else (None, None)
            if not min_next: min_next = min_each
            if not max_next: max_next = max_each
            # identify pairs
            if int(max_each) < int(min_next) or int(min_each) > int(max_next):
                get_minimum_maximum = []
                rln = 0
                if r11 == None and r22 == None:
                    get_minimum_maximum = [int(r1), int(r2)]
                    rln = abs(int(r1) - int(r2))
                else:
                    get_minimum_maximum = [int(r1), int(r2), int(r11), int(r22)]
                    rln = abs(int(r1) - int(r2)) + abs(int(r11) - int(r22))
                minr, maxr = min(get_minimum_maximum), max(get_minimum_maximum)
                junction_looper(minr, maxr, rln, junctions, key, key2, cp_key2)
            else:
                rln = abs(int(r1) - int(r2))
                junction_looper(r1, r2, rln, junctions, key, key2, cp_key2)
                continue

with open(os.path.join(output_dir, f"{id}.readHashc.pkl"), "wb") as f:
    pickle.dump(dict(readHashc), f)