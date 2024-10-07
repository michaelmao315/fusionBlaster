import pickle
import sys
import os
import glob

fastq = sys.argv[1]
output_dir = sys.argv[2]

readHashc_combined = {}

# Get all readHashc files for this fastq
readHashc_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*_*.readHashc.pkl"))

for file_path in readHashc_files:
    with open(file_path, "rb") as f:
        try:
            readHashc_part = pickle.load(f)
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
            continue

        for key in readHashc_part:
            if key not in readHashc_combined:
                readHashc_combined[key] = {}
            for key2 in readHashc_part[key]:
                if key2 not in readHashc_combined[key]:
                    readHashc_combined[key][key2] = {}
                for key3 in readHashc_part[key][key2]:
                    readHashc_combined[key][key2][key3] = readHashc_combined[key][key2].get(key3, 0) + readHashc_part[key][key2][key3]

output_file = os.path.join(output_dir, f"{fastq}.readHashc.pkl")
print(f"Writing combined data to: {output_file}")  # Debug print
with open(output_file, "wb") as f:
    pickle.dump(readHashc_combined, f)

print("Combination of readHashc complete.")

for file_path in glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*_*.readHashc.pkl")):
    try:
        os.remove(file_path)
    except Exception as e:
        print(f"Error removing {file_path}: {str(e)}")
print("Cleanup complete.")