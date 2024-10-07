import os
import sys

def find_cut_position(lines, target_lines):
    current_read_name = lines[-1].split('\t')[0]
    cut_position = len(lines) - 1

    while cut_position >= 0:
        read_name = lines[cut_position].split('\t')[0]
        if read_name != current_read_name:
            break
        cut_position -= 1

    return cut_position + 1

def split_file(input_file, output_prefix, target_lines):
    file_count = 0
    chunk = []

    with open(input_file, 'r') as f:
        for line in f:
            chunk.append(line)

            if len(chunk) >= target_lines:
                cut_position = find_cut_position(chunk, target_lines)
                file_count += 1
                output_file = f"{output_prefix}_{file_count:04d}.tsv"

                with open(output_file, 'w') as out_f:
                    out_f.writelines(chunk[:cut_position])

                chunk = chunk[cut_position:]

        if chunk:
            file_count += 1
            output_file = f"{output_prefix}_{file_count:04d}.tsv"

            with open(output_file, 'w') as out_f:
                out_f.writelines(chunk)

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    target_lines = int(sys.argv[3])

    split_file(input_file, output_prefix, target_lines)