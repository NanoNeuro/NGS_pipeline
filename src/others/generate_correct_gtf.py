import argparse

def add_transcript_id_version(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # Skip comment lines
                outfile.write(line)
                continue
            
            columns = line.strip().split('\t')
            
            column_text = columns[-1]
            columns_dict = {i.split(' "')[0]: i.split(' "')[1][:-1] for i in column_text.split("; ")}
            
            
            if ("transcript_id" in columns_dict.keys()) & ("transcript_version" in columns_dict.keys()):
            	transcript_id = columns_dict['transcript_id']
            	transcript_version = columns_dict['transcript_version']
            	
            	columns[-1] += f' transcript_id_version "{transcript_id}.{transcript_version}";'

            outfile.write('\t'.join(columns) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add 'transcript_id_version' attribute to GTF file.")
    parser.add_argument('--input', required=True, help='Input GTF file path')
    parser.add_argument('--output', required=True, help='Output GTF file path')

    args = parser.parse_args()
    add_transcript_id_version(args.input, args.output)

    print(f"New GTF file created: {args.output}")
