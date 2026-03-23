using FASTX

function cal_N_ins(Indelprofile::String, Site::Int64)
    N_ins = Dict{String,Vector{Int64}}()
    bases = ["A", "C", "G", "T"]
    #bases = [x for x in bases if x != refbase]
    fills = (string(Site + 10) * "I1:") .* bases
    #push!(fills, string(Site + 9) * "I1:" * refbase)
    for fill in fills
        N_ins[fill] = Vector{Int64}()
    end
    for line in eachline(Indelprofile)
        if startswith(line, "reference")
            continue
        else
            _, _, count, insertions, deletions, _ = split(line, "\t")
            ct = parse(Int64, count)
            for tp in fills
                if occursin(tp, insertions)
                    push!(N_ins[tp], ct)
                end
            end
        end
    end
    Dict(zip(keys(N_ins), map(i -> (sum(i)), values(N_ins))))
end

function read_indel_file(indelFile::String)
    indels = Dict{String,Tuple{Int64,Int64}}()
    for line in eachline(indelFile)
        if startswith(line, "gRNAID")
            continue
        else
            libraryID, total, indel = split(line, "\t")
            libraryID = string(libraryID)
            indels[libraryID] = (parse(Int64, total), parse(Int64, indel))
        end
    end
    indels
end


function fin(input::String, indelFile::String, Site::Int64, reference::String, output::String)
    refs = open(FASTA.Reader, reference, index="$reference.fai")
    all_lib = collect(keys(refs.index.names))
    all_indels = read_indel_file(indelFile)
    open("$output/N_$(Site)_ins_fullgRNA_30mer.txt", "w") do io
        for i in all_lib
            N_seq = FASTX.sequence(refs[i])[110:139]
            N_base = string(N_seq[Site+4])
            if isfile("$input/$i.indelProfile.tsv") && stat("$input/$i.indelProfile.tsv").size > 0
                N_ins = cal_N_ins("$input/$i.indelProfile.tsv", Site)
                totals, N_indel = all_indels[i]
                for (k, v) in N_ins
                    k = replace(k,r"\d+I1:"=>"")
                    println(io, i, "\t", N_seq, "\t", k, "\t", v, "\t", N_indel, "\t",totals, "\t", N_base)
                end
            end
        end
    end
end

function main()
    Threads.@threads for tp in ["P1305-1","P1305-2","P1307-1","P1337-2","P1339-1","P1339-2"]
       fin("indel_profile/HEK_on_$(tp)","indelCounts/HEK_on_$(tp).indel_efficiency.tsv",17,"refs/ref_on_144.fa","$(tp)_1bpIns")
    end
end

main()
