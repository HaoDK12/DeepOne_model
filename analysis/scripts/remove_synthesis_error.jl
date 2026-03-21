using ArgMacros
using FASTX

function ReadIn(fs::String)
    seq_lst = String[]
    for line in eachline(fs)
        push!(seq_lst, line)
    end
    seq_lst
end


function removeBAC(ref::String, ctrl::Vector{String}, case::Vector{String}, start::Int64, fin::Int64)
    filter(x -> (!((x in ctrl) && x != ref[start:fin])), case)
end

function fin()
    @inlinearguments begin
                @argumentrequired String spCas9 "-s" "--spCas9"
                @argumentrequired String WT "-w" "--WT"
                @argumentrequired String reference "-r" "--reference"
                @argumentrequired String outPath "-o" "--outPath"
                @argumentrequired Int64 start "-st" "--start"
                @argumentrequired Int64 fin "-f" "--fin"
        end
    
    
    refs = open(FASTA.Reader,reference)
    if isdir(outPath) == false
        mkdir(outPath)
    end
    for record in refs
        refname = identifier(record)
        refseq = convert(String,sequence(record))
        if isfile("$(spCas9)/$(refname).out") && isfile("$(WT)/$(refname).out")
            case = ReadIn("$(spCas9)/$(refname).out")
            ctrl = ReadIn("$(WT)/$(refname).out")
            open("$outPath/$(refname).filter.out", "w") do IO
                for line in removeBAC(refseq, ctrl, case, start, fin)
                    write(IO, line, "\n")
                end
            end
        else
            println("$refname have been removed!")
        end
    end
end

fin()
