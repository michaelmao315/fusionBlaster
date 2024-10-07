import pickle
import sys
import os
import glob

fastq = sys.argv[1]
output_dir = sys.argv[2]

readHash_combined = {}

# Get all readHash files for this fastq
readHash_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*_*.readHash.pkl"))

for file_path in readHash_files:
    with open(file_path, "rb") as f:
        try:
            readHash_part = pickle.load(f)
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
            continue

        for key in readHash_part:
            if key not in readHash_combined:
                readHash_combined[key] = {}
            for key2 in readHash_part[key]:
                if key2 not in readHash_combined[key]:
                    readHash_combined[key][key2] = {}
                for key3 in readHash_part[key][key2]:
                    readHash_combined[key][key2][key3] = readHash_combined[key][key2].get(key3, 0) + readHash_part[key][key2][key3]

output_file = os.path.join(output_dir, f"{fastq}.readHash.pkl")
print(f"Writing combined data to: {output_file}")  # Debug print
with open(output_file, "wb") as f:
    pickle.dump(readHash_combined, f)

print("Combination of readHash complete.")

for file_path in glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*_*.readHash.pkl")):
    try:
        os.remove(file_path)
    except Exception as e:
        print(f"Error removing {file_path}: {str(e)}")

print("Cleanup complete.")