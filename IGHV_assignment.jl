## Script to analyse and identify IGHV sequences ##

# Alaine Athenaïs Marsden - 2024: alainem@nicd.ac.za #

################################################################################################################

## Require Pre-Processing Steps ##


#= This script requires several inputs:

1. preprocessed participant sequences - qc'ed, trimmed, dereplicated
2. reference IGHV region sequences: incl -> v-region, RSS, LPART1 and LPART2
3. Tables generated from local BLAST with custom v-region db

script assumes a directory full of participant specific directories, aptly named.

=#

################################################################################################################

## Necessary Libraries ##

using BioAlignments, CSV, DataFrames, FASTX, StatsBase

################################################################################################################

## Defined Functions and Structures ##

# function to compile sequences from fastas
function compile_seqs(file)
    seqs = FASTA.Reader(open(file, "r"))
    seq_list = []
    for i in seqs
        push!(seq_list, i)
    end
    return seq_list
end

# function to assign IGHV alleles through parsing blast tab
function assign_alls(participant_seqs, blast_tab)
    assigned_seqs = []
    for i in participant_seqs
        seq_id = identifier(i)
        seq = sequence(i)
        seq_tab = filter(:qseqid => ==(seq_id), blast_tab)
        #find best match with lowest evalue
        min_eval = minimum(seq_tab.evalue)
        filter!(:evalue => ==(min_eval), seq_tab)
        #some sequences may need further processing via pident 
        max_pident = maximum(seq_tab.pident)
        filter!(:pident => ==(max_pident), seq_tab)
        #get best match
        id_row = first(seq_tab,1)
        all_match = id_row.sseqid
        if id_row.mismatch > 0
            all_match = join([split(all_match, "*")[1], "novel"], "*")
        end
        seq_id = join([seq_id, all_match], "_")
        new_record = FASTA.Record(seq_id, seq)
        push!(assigned_seqs, new_record)
    end
    return assigned_seqs
end

# genotype data_structure

struct genotype
    cap
    gene_dict
end

# lists to fasta sequences

function seq_gen(seq_list, id_list)
    seq_records = []
    for i in 1:eachindex(seq_list)
        seq = seq_list[i]
        id = id_list[i]
        seq_rec = FASTA.Record(id, seq)
        push!(seq_records, seq_rec)
    end
    return seq_records
end

# fasta file generator
function fasta_gen(dir, rec_list)
    x = open(FASTA.Writer, dir)
    for i in rec_list
        write(x, i)
    end
    close(x)
end

################################################################################################################

## Load In Relevant Data, insert locations instead of 'X' ##

#reference data

v_region_refs = compile_seqs(x) 
lpart1_refs = compile_seqs(x)
lpart2_refs = compile_seqs(x)
rss_refs = compile_seqs(x)

#directory with participant directories

working_directory = x
cd(working_directory)

# get participant folders - adjust as needed, written for CAP samples

cap_dirs = []
for i in readdir()
    if startswith(i, "CAP")
        push!(cap_dirs,i)
    end
end

################################################################################################################

## Begin iterating through participant directories and determining participant genotypes and compiling alleles for each gene ##
# replace x with approriate file_name conventions

cap_genotypes = []
alls_per_gene = Dict()

for i in cap_dirs
    cap_tab = CSV.read(open(x), DataFrame) #blast output tab
    cap_seqs = compile_seqs(x) #processed participant seqs
    assigned_cap_seqs = assign_alls(cap_seqs, cap_tab)
    # create genotype
    gene_dict = Dict()
    for seq in assigned_cap_seqs
        seq_id = identifier(i)
        gene = split(split(seq_id, "_")[2], "*")[1] #these numbers will vary based on what other info is in sequence id: size, readnumber etc
        if gene ∉ keys(gene_dict)
            gene_dict[gene] = [seq]
        else
            push!(gene_dict[gene], seq)
        end
    end
    cap_geno = genotype(i, gene_dict)
    push!(cap_genotypes, cap_geno)
    # begin compiling sequences per gene
    for seq in assigned_cap_seqs
        seq_id = identifier(i)
        gene = split(split(seq_id, "_")[2], "*")[1] #these numbers will vary based on what other info is in sequence id: size, readnumber etc
        if gene ∉ keys(alls_per_gene)
            alls_per_gene[gene] = [seq]
        else
            push!(alls_per_gene[gene], seq)
        end
    end
end

################################################################################################################

## Get unique alleles ##

gene_list = keys(alls_per_gene)
unique_seqs = []
unique_ids = []
for gene in gene_list
    gene_alls = alls_per_gene[gene]

    for a in gene_alls
        if sequence(a) ∉ unique_seqs
            push!(unique_seqs, sequence(a))
            new_id = split(identifier(a), "_")[2]*"_"*split(identifier(a), "_")[3]
            push!(unique_ids,new_id)
        end
    end
end

unique_alleles = seq_gen(unique_seqs, unique_ids)
fasta_gen(X, unique_alleles) #replace x with output directory

# determine duplications in IGHV germline repertoire

cap_cnv = Dict()

for gene in gene_list
    for cap in cap_genotypes
        cnv = "none"
        cap_alls = cap.gene_dict[gene]
        if length(cap_alls) > 1
            cnv = "dup"
        end
        cap_cnv[cap] = cnv
    end
end

