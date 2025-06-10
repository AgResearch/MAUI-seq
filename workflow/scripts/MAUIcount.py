#!/usr/bin/env python3
import os
import subprocess
import datetime
import yaml
import argparse

mauiseq_description = """
MAUIcount: Analyze amplicon sequences with UMI-based error correction

This script reads amplicon sequences from a set of fastq files.
The first UMI_len bases of each read are a random tag (UMI or Unique Molecular Identifier).
The script keeps track of how many times each UMI is used with each unique sequence.
For the set of samples, outputs are files with a list of fasta sequences in descending rank
order of abundance, and corresponding tables with the counts of each sequence in each sample.

Three output sets are produced by default:

"accepted_sequences": the main MAUI-seq output, using secondary sequences to filter out errors

"all_primary_sequences": the same UMI counts before error filtering

"read_sequences": 'conventional' analysis of the sequences, ignoring UMIs

An additional output set can be produced by setting output_types = 4. This has the same 
filtering as accepted_sequences, but applied on a per-sample basis, rather than on totals 
across all samples. If allele frequencies vary greatly across samples, this would be 
preferable in principle, but it is not reliable unless read counts are very high. Otherwise,
it can lead to sequences being stochastically deleted from some samples but not others.

In all cases, the outputs are truncated to discard rare sequences that would have
frequencies below add_limit in the overall set of samples. 

Additional output files are summary.txt (various data about the run) and UMI_stats.tab 
(information on the distribution of read numbers per UMI that may be useful for optimising
the protocol).

A file amplicon_config.yaml specifies the files to be analysed and gene-specific parameters.

If there is a file fastq_file_list.txt in the same folder as the data, only the files listed 
in this file will be included in the analysis. If this file is not present, it will be  
created with a list of all files that have the extension .fastq.

"""
epilog_text ="""
Written by Peter Young. Version 1.01 on 20 January 2020.
Refactored by Ben Perry. Version 2.0 on 10 June 2025.
"""

start_time = datetime.datetime.now()
script_file = os.path.basename(__file__)

#Parse command line arguments
def parse_arguments():
    """
    Parse command line arguments for MAUIcount script.
    """
    
    parser = argparse.ArgumentParser(
        description=mauiseq_description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog_text
    )
    
    parser.add_argument(
        "--working-folder",
        type=str,
        default="./",
        help="Path to folder containing FASTQ files."
    )
    
    parser.add_argument(
        "--amplicon-config",
        type=str,
        default="resources/amplicon_config.yaml",
        help="Path to amplicon configuration file."
    )
    
    parser.add_argument(
        "--fastq-file-list",
        type=str,
        required=False,
        default=None,
        help="Path to a file listing FASTQ files to process. If not provided, all .fastq files in the working folder will be used."
    )

    parser.add_argument(
        "--read-diff",
        type=int,
        default=2,
        help="Count UMI only if most abundant sequence has at least read_diff more reads than the next."
    )
    
    parser.add_argument(
        "--reject-threshold",
        type=float,
        default=0.7,
        help="Reject sequences that occur as second sequences with UMIs at least reject_threshold times as often as they occur as the primary sequence."
    )
    
    parser.add_argument(
        "--add-limit",
        type=float,
        default=0.001,
        help="Sequences are included, in rank order, until the next would add a fraction less than add_limit (i.e. this discards sequences with an overall relative abundance less than add_limit)."
    )
    
    parser.add_argument(
        "--output-types",
        type=int,
        default=3,
        choices=[3, 4],
        help="If set to 4, additional output is produced based on filtering separately for each sample."
    )
    
    parser.add_argument(
        "--output-read-counts",
        action="store_true",
        help="Output the total number of reads contributing to the UMI primary and secondary counts."
    )
    
    parser.add_argument(
        "--gene",
        type=str,
        required=True,
        choices=["recA", "rpoB", "nodA", "nodD", "viciae_nodD"],
        default="recA",
        help="Gene amplicon to analyze."
    )
    
    return parser.parse_args()

#Parse command line arguments
args = parse_arguments()

#Functions
def load_amplicon_parameters(config_file="resources/amplicon_config.yaml", gene_name="recA"):
    """
    Load amplicon-specific parameters from resources/amplicon_parameters.yaml configuration file.
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config['genes'][gene_name]

def find_match(line,dic):
    """
    Split a sequence line into UMI and sequence (removing the primers);
    then add to a dictionary of all the different UMIs with counts of each of their sequences
    Uses global constants UMI_len, f_primer_len, r_primer_len
    Calls increment()
    """
    global UMI_len, f_primer_len, r_primer_len

    UMI = line[0:UMI_len]
    sequence = line[(UMI_len + f_primer_len):(len(line) - r_primer_len)]
    if UMI in dic:
        increment(dic[UMI],sequence,1)
    else:
        dic[UMI] = {sequence:1}


def increment(dic,key,count):
    """
    Increments the total for key in dic by count, or creates key:count if not existing
    """
    if key in dic:
        dic[key] += count
    else:
        dic[key] = count    



def main():
    """
    Main function to run the MAUIcount script.
    """
    # amplicon_config.yaml parameters
    global UMI_len, f_primer_len, r_primer_len, total_len
    # Command line arguments
    global read_diff, reject_threshold, add_limit, output_types, output_read_counts, working_folder

    #Initialise lists and dictionaries
    sample_tuples = []
    all_samples_table = {}
    total_counts = {}
    total_clean_by_sample_counts = {}
    clean_all_samples_table = {}
    total_sec_seq_counts = {} 
    total_cleaned_sec_counts = {}
    all_samples_reads = {}
    total_reads = {}
    raw_read_count = 0
    total_pri_reads = {} #NEW
    total_sec_reads = {} #NEW
    UMI_stat_names = ["all sequences per UMI", "primary sequences per UMI", 
    "secondary sequences per UMI", "distinct sequences per UMI", "pri - sec difference",
    "UMI profile"]
    UMI_stats = [{},{},{},{},{},{}] #NEW
    
    # Map command line arguments to global variables
    read_diff = args.read_diff
    reject_threshold = args.reject_threshold
    add_limit = args.add_limit
    output_types = args.output_types
    output_read_counts = args.output_read_counts
    working_folder = args.working_folder
    
    # Load parameters from amplicon config
    try:
        amplicon_config = load_amplicon_parameters(args.amplicon_config, args.gene)
    except (FileNotFoundError, KeyError) as e:
        print(f"Error: Could not load amplicon parameters from {args.config}: {e}")
        return
    
    # Extract amplicon-specific parameters from config
    UMI_len = amplicon_config['UMI_len']
    f_primer_len = amplicon_config['f_primer_len']
    r_primer_len = amplicon_config['r_primer_len']
    total_len = amplicon_config['total_len']


    # Set working folder
    if not os.path.isdir(working_folder):
        print("Error: The specified working folder does not exist.")
        return


    #Get a list of the samples to process
    #Use the list file if it already exists, or else create one
    fastq_path = os.path.join(working_folder, "*.fastq")
    if not os.path.isfile(fastq_path):
        fastq_list_path = os.path.join(working_folder, "fastq_file_list.txt")
        subprocess.call("ls " + fastq_path + " > " + fastq_list_path, shell=True)

    #The sample ID will be the fasta filename up to the first "."
    with open(fastq_list_path) as file_list:
        for fastq_filename in file_list:

            if "/" in fastq_filename:
                fastq_filename = fastq_filename[fastq_filename.rindex("/")+1:]
            
            fastq_filepath = os.path.join(working_folder, fastq_filename.rstrip("\n"))
            sample_ID = fastq_filename[0:fastq_filename.find(".")]
            sample_tuples.append((sample_ID,fastq_filepath))

                
    #Process the read data, one sample at a time

    for (sample_ID,fastq_filepath) in sample_tuples:                        
        UMI_dict = {}
        read_dict = {} #keys will be sequences, items will be number of reads for each sequence in the sample

        #Read in a fastq file one sequence record at a time (4 lines) and process the DNA sequence (2nd line) with find_match
        #Create UMI_dict with structure {UMI:{sequence:count,...},...}
        with open(fastq_filepath) as fastq_file:
            ctr = 0
            record = []
            for next_line in fastq_file:
                record.append(next_line.rstrip("\n"))
                ctr += 1
                if ctr == 4:
                    raw_read_count +=1
                    if len(record[1]) == total_len: #Only process sequences that are the correct length
                        find_match(record[1],UMI_dict)
                        sequence = record[1][(UMI_len + f_primer_len):(len(record[1]) - r_primer_len)]
                        increment(read_dict, sequence, 1) 
                        increment(total_reads, sequence, 1)
                    record = []
                    ctr = 0

        #Extract sequence count data from UMI_dict
        sample_counts = {}  #keys will be sequences, items will be number of UMIs for each sequence in the sample
        sec_seq_counts = {} #keys will be sequences, items will be how often they are second sequence in a UMI           
        clean_sample_counts = {} #only includes sequences that are below the threshold for second sequence count
            
        for UMI, matches in sorted(UMI_dict.items(), key=lambda item: sum(item[1].values()), reverse=True):
    
            sorted_list = sorted(UMI_dict[UMI].items(), key=lambda item: item[1], reverse=True)
            sequence = sorted_list[0] #Choose the most abundant sequence for each UMI
            if len(sorted_list) > 1: #There is a second sequence with this UMI
                sec_seq = sorted_list[1]
            else:
                sec_seq = ("null",0)
                
            if sequence[1] - sec_seq[1] >=read_diff: 
                #Only include sequences that have at least read_diff more reads than the sec_seq                  
                increment(sample_counts, sequence[0], 1)
                increment(total_counts, sequence[0], 1)
                increment(total_pri_reads, sequence[0], sequence[1]) #NEW
                
                #Count second sequences across all UMIs
                if sec_seq[0] != "null":
                    increment(sec_seq_counts, sec_seq[0], 1)
                    increment(total_sec_seq_counts, sec_seq[0], 1)
                    increment(total_sec_reads, sec_seq[0], sec_seq[1]) #NEW
                        
            #NEW: Collect reads-per-UMI distribution
            increment(UMI_stats[0], sum(UMI_dict[UMI].values()), 1)
            increment(UMI_stats[1], sequence[1], 1)
            increment(UMI_stats[2], sec_seq[1], 1)
            increment(UMI_stats[3], len(UMI_dict[UMI]), 1)
            increment(UMI_stats[4], sequence[1] - sec_seq[1], 1)
            count_string = ""
            for sequence, count in sorted_list:
                count_string += str(count) + ","
            increment(UMI_stats[5], count_string, 1)
            
                            
        for seq in sample_counts:
            seq_count = sample_counts[seq]
            if seq in sec_seq_counts: 
                sec_count = sec_seq_counts[seq]
            else: sec_count = 0
            
            if sec_count < seq_count * reject_threshold: 
                clean_sample_counts[seq] = seq_count
                increment(total_clean_by_sample_counts, seq, seq_count)
                increment(total_cleaned_sec_counts, seq, sec_count)

        all_samples_table[sample_ID] = sample_counts
        clean_all_samples_table[sample_ID] = clean_sample_counts    
        
        all_samples_reads[sample_ID] = read_dict #Record read data for this sample
        
        

        
    #Prepare to output the results

    output_folder = os.path.join(working_folder , "MAUIcount_output") 
    os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
    

    fas_file = {}
    table_file = {}                
    file_names = ["all_primary_sequences","accepted_sequences","read_sequences","accepted_by_sample_sequences"]
    for i in range(output_types):
        fas_file[i] = open(os.path.join(output_folder, (file_names[i] + ".fas")), "w")
        table_file[i] = open(os.path.join(output_folder, (file_names[i] + ".tab")), "w")

    rank = 0
    ranked_sequence_list = []
    sequence_list = [[],[],[],[]]
    cumultotal = [0,0,0,0]
    total_accepted_seqs = 0
    total_accepted_counts = 0
    reported_UMI_seqs = 0
    reported_accepted_seqs = 0
    reported_read_seqs = 0

    
    #Write sets of sequences in fasta format
    #Headers have raw sequence rank, total counts, total secondary counts  

    for sequence, seq_count in sorted(total_counts.items(), key=lambda item: item[1], reverse=True):
        
        if sequence in total_sec_seq_counts:
            sec_seq_count = total_sec_seq_counts[sequence]
        else: sec_seq_count = 0
    
        rank += 1
        ranked_sequence_list.append(sequence)
        
        if seq_count > cumultotal[0]*add_limit:
            cumultotal[0] += seq_count
            sequence_list[0].append(sequence)
            reported_UMI_seqs +=1
            
            fas_file[0].write(">seq_%d_%d_%d\n%s\n" % (rank, seq_count, sec_seq_count, sequence))
            
        if sec_seq_count < seq_count*reject_threshold:
            total_accepted_seqs += 1
            total_accepted_counts += seq_count 
            if seq_count > cumultotal[1]*add_limit:
                cumultotal[1] += seq_count
                sequence_list[1].append(sequence)
                reported_accepted_seqs +=1
                
                fas_file[1].write(">seq_%d_%d_%d\n%s\n" % (rank, seq_count, sec_seq_count, sequence))

    #if required, output a similar file filtered by treating each sample separately        
    if output_types == 4: 
        for sequence, seq_count in sorted(total_clean_by_sample_counts.items(), key=lambda item: item[1], reverse=True):
            if seq_count > cumultotal[3]*add_limit:   #omit sequences with few counts
                cumultotal[3] += seq_count
            
                if sequence in total_cleaned_sec_counts:
                    sec_seq_count = total_cleaned_sec_counts[sequence]
                else: sec_seq_count = 0
        
                rank = ranked_sequence_list.index(sequence) + 1
                sequence_list[3].append(sequence)
            
                fas_file[3].write(">seq_%d_%d_%d\n%s\n" % (rank, seq_count, sec_seq_count, sequence))
            
    #Raw reads are ranked and named differently (seqr) because not all are in the UMI list 
    #Headers have raw sequence rank, total counts                          
    rank_r = 0
    ranked_sequence_list_r = []

    for sequence, seq_count in sorted(total_reads.items(), key=lambda item: item[1], reverse=True):
        rank_r +=1
        ranked_sequence_list_r.append(sequence)    
        if seq_count > cumultotal[2]*add_limit:
            cumultotal[2] += seq_count
            sequence_list[2].append(sequence)
            reported_read_seqs +=1
            fas_file[2].write(">seqr_%d_%d\n%s\n" % (rank_r, seq_count, sequence))

        
    #Write tables of counts for each sequence in each sample
            
    #Write sequence ranks as column headers
    for sequence in sequence_list[0]:
        rank = ranked_sequence_list.index(sequence) +1
        table_file[0].write("\tseq_%d" % (rank))
    table_file[0].write("\n")    
    for sequence in sequence_list[1]:
        rank = ranked_sequence_list.index(sequence) +1
        table_file[1].write("\tseq_%d" % (rank))
    table_file[1].write("\n")

    #Write a row of counts for each sample
    for sample_ID in sorted(all_samples_table):
        for i in range(2):
            table_file[i].write(sample_ID+"\t")        
        for sequence in ranked_sequence_list:
            if sequence in all_samples_table[sample_ID]:
                seq_count = all_samples_table[sample_ID][sequence]
            else:
                seq_count = 0
            for i in range(2):
                if sequence in sequence_list[i]:    
                    table_file[i].write(str(seq_count)+"\t")                    
        for i in range(2):
            table_file[i].write("\n")        

    #Write total counts for each sequence under corresponding column
    for i in range(2):
        table_file[i].write("total")
        for sequence in sequence_list[i]:
            table_file[i].write("\t"+str(total_counts[sequence]))
        table_file[i].write("\nseconds")
        for sequence in sequence_list[i]:
            if sequence in total_sec_seq_counts:
                sec_seq_count = total_sec_seq_counts[sequence]
            else: sec_seq_count = 0    
            table_file[i].write("\t"+str(sec_seq_count))

        #if required, output pri/sec read counts
        if output_read_counts:        
            table_file[i].write("\n\npri reads")
            for sequence in sequence_list[i]:
                table_file[i].write("\t"+str(total_pri_reads[sequence]))
            table_file[i].write("\nsec reads")
            for sequence in sequence_list[i]:
                if sequence in total_sec_reads:
                    sec_reads = total_sec_reads[sequence]
                else: sec_reads = 0    
                table_file[i].write("\t"+str(sec_reads))
            
    #if required, output a similar table filtering each sample separately        
    if output_types == 4:
        for sequence in sequence_list[3]:
            rank = ranked_sequence_list.index(sequence) +1
            table_file[3].write("\tseq_%d" % (rank))
        table_file[3].write("\n")

        for sample_ID in sorted(clean_all_samples_table):
            table_file[3].write(sample_ID+"\t")        
            for sequence in sequence_list[3]:
                if sequence in clean_all_samples_table[sample_ID]:
                    seq_count = clean_all_samples_table[sample_ID][sequence]
                else:
                    seq_count = 0
                table_file[3].write(str(seq_count)+"\t")                            
            table_file[3].write("\n")

        table_file[3].write("total")
        for sequence in sequence_list[3]:
            table_file[3].write("\t"+str(total_clean_by_sample_counts[sequence]))
        table_file[3].write("\nseconds")
        for sequence in sequence_list[3]:
            if sequence in total_cleaned_sec_counts:
                sec_seq_count = total_cleaned_sec_counts[sequence]
            else: sec_seq_count = 0    
            table_file[3].write("\t"+str(sec_seq_count))

    #Write a similar table based on reads

    #Write sequence ranks as column headers
    for sequence in sequence_list[2]:
        rank_r = ranked_sequence_list_r.index(sequence) +1
        table_file[2].write("\tseqr_%d" % (rank_r)) 
    table_file[2].write("\n")

    #Write a row of counts for each sample
    for sample_ID in sorted(all_samples_reads):
        table_file[2].write(sample_ID+"\t")    
        for sequence in sequence_list[2]:
            if sequence in all_samples_reads[sample_ID]:
                seq_count = all_samples_reads[sample_ID][sequence]
            else:
                seq_count = 0
            table_file[2].write(str(seq_count)+"\t")                    
        table_file[2].write("\n")

    #Write total counts for each sequence under corresponding column
    table_file[2].write("total")
    for sequence in sequence_list[2]:
        table_file[2].write("\t"+str(total_reads[sequence]))
        
    #Close all output files
    for i in range(output_types):
        fas_file[i].close()
        table_file[i].close()
        
    #NEW: Write a table with UMI statistics
    with open(os.path.join(output_folder, "UMI_stats.tab"), "w") as stats:
        for i in range(len(UMI_stats)):
            stats.write("\n" + UMI_stat_names[i] + "\n")
            for key,value in sorted(UMI_stats[i].items()):
                stats.write(str(key) + "\t" + str(value) + "\n")

    #Write a text file with a summary of parameters and various counts
    with open(os.path.join(output_folder, "summary.txt"), "w") as summary:
        summary.write(script_file +"\n")
        summary.write(working_folder + "\n\n")
        summary.write(str(read_diff) + "\tread_diff\n")
        summary.write(str(reject_threshold) + "\treject_threshold\n")
        summary.write(str(add_limit) + "\tadd_limit\n\n")
        summary.write(str(raw_read_count) + "\tRaw read count\n")
        summary.write(str(sum(total_reads.values())) + "\tTotal reads used\n")
        summary.write(str(len(total_reads)) + "\tTotal unique sequences\n")
        summary.write("\nBefore truncation:\n")
        summary.write(str(sum(total_counts.values())) + "\tTotal UMI counts\n")
        summary.write(str(total_accepted_counts) + "\tTotal accepted counts\n")
        summary.write(str(len(total_counts)) + "\tTotal UMI sequences\n")
        summary.write(str(total_accepted_seqs) + "\tTotal accepted sequences\n")
        summary.write("\nAfter truncating to threshold:\n")
        summary.write(str(reported_UMI_seqs) + "\tReported UMI sequences\n")
        summary.write(str(reported_accepted_seqs) + "\tReported accepted sequences\n\n")
        summary.write(str(reported_read_seqs) + "\tReported read sequences\n")
    
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        summary.write("\n")
        summary.write(str(start_time) + "\tStart time\n")
        summary.write(str(end_time) + "\tEnd time\n")
        summary.write(str(duration) + "\tDuration")
        
    print("MAUIcount completed in " + str(duration.total_seconds()) + " seconds")


if __name__ == "__main__":
    main()