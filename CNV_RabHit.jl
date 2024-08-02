## Script to analyse and identify IGHV sequences ##

# Alaine Athenaïs Marsden - 2024: alainem@nicd.ac.za #

################################################################################################################

#= This script requires the tabular output from Rabhit and outputs a table calculating fold change in expression
with duplication =#

################################################################################################################

## Necessary Libraries ##

using CSV, DataFrames, StatsPlots, Statistics

################################################################################################################

halpo = CSV.read(open(X), DataFrame, delim = "\t") #input rabhit table where x

# Tidy up table 

cols = names(haplo)

filter!(x->x≠"Column16",cols)
select!(haplo, Not("subject"))
rename!(haplo, cols)

replace!(haplo.counts1, "NA" => "0,0")
replace!(haplo.counts2, "NA" => "0,0")
replace!(haplo.counts3, "NA" => "0,0")
replace!(haplo.counts4, "NA" => "0,0")
transform!(haplo, :counts1 => ByRow(x -> split(x, ',')) => [:ch1c1, :ch2c1])
transform!(haplo, :counts2 => ByRow(x -> split(x, ',')) => [:ch1c2, :ch2c2])
transform!(haplo, :counts3 => ByRow(x -> split(x, ',')) => [:ch1c3, :ch2c3])
transform!(haplo, :counts4 => ByRow(x -> split(x, ',')) => [:ch1c4, :ch2c4])
haplo.ch1c1 = parse.(Int64, haplo.ch1c1)
haplo.ch1c2 = parse.(Int64, haplo.ch1c2)
haplo.ch1c3 = parse.(Int64, haplo.ch1c3)
haplo.ch1c4 = parse.(Int64, haplo.ch1c4)
haplo.ch2c1 = parse.(Int64, haplo.ch2c1)
haplo.ch2c2 = parse.(Int64, haplo.ch2c2)
haplo.ch2c3 = parse.(Int64, haplo.ch2c3)
haplo.ch2c4 = parse.(Int64, haplo.ch2c4)

################################################################################################################

# calculate proportion for each gene and determine CNV type

proportion_count = []
cnv = []
by_gene = []

for i in eachrow(haplo)

    d1 = false
    del = false
    cap_total = cap_tot[i.subject]
    g_mean = gene_means[i.gene]
    g_range = gene_ranges[i.gene]
    g_std = gene_stds[i.gene]
    

    if i["IGHJ6_02"] == "Unk"
        continue
    end

    push!(by_gene, i.gene)

    if length(split(i["IGHJ6_02"], ",")) > 1 
        d1 = true
    end

    if i["IGHJ6_02"] == "Del" 
        del = true
    end

    if del
        push!(cnv, "del")
    elseif d1 
        push!(cnv, "dup")
    else
        push!(cnv, "none")
    end

    count = sum([i.ch1c1, i.ch1c2, i.ch1c3, i.ch1c4])/cap_total
    push!(proportion_count, count)
end

for i in eachrow(haplo)

    if i.gene == "IGHV3-64D" || i.gene == "IGHV2-70D" || i.gene == "IGHV4-28" || i.gene == "IGHV1-45"
        continue
    end

    d2 = false
    del = false
    

    if i["IGHJ6_02"] == "Unk"
        continue
    end

    push!(by_gene, i.gene)

    if length(split(i["IGHJ6_02"], ",")) > 1 
        d2 = true
    end

    if i["IGHD2-8_02"] == "Del" 
        del = true
    end

    if del
        push!(cnv, "del")
    elseif d2
        push!(cnv, "dup")
    else
        push!(cnv, "none")
    end

    count = sum([i.ch2c1, i.ch2c2, i.ch2c3, i.ch2c4])/cap_total
    push!(proportion_count, count)
end

df = DataFrame(cnv=cnv, gene=by_gene, proportion=proportion_count)

################################################################################################################

# Determine expression changes

CSV.write(X, by_cnv) #set output directory and filename where X


dup_genes = []
ave_dup = []
ave_none = []
fold_change = []

for k in genes
    gene_tab = filter(:gene => isequal(k), df)
    if "dup" ∈ gene_tab.cnv
        gene_tab = filter(:cnv => !isequal("del"), gene_tab)
        none_tab = filter(:cnv => !isequal("none"), gene_tab)
        dup_tab = filter(:cnv => !isequal("dup"), gene_tab)

        ave_d = mean(log.(dup_tab.proportion))
        ave_n = mean(log.(none_tab.proportion))
        fold = (ave_d/ave_n)-1

        push!(ave_dup, ave_d)
        push!(ave_none, ave_n)
        push!(fold_change, fold)
        push!(dup_genes, k)
    end
end

fold_tab = DataFrame()
fold_tab[!, "gene"] = dup_genes
fold_tab[!,"fold"] = fold_change

CSV.write(X, fold_tab) #set output directory and filename where X