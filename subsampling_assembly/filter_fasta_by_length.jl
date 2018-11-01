#!/usr/bin/env julia

using Pkg
using FastaIO
using ArgParse
using Logging

# using FastaIO 
#try
#    Pkg.installed("FastaIO")
#    using FastaIO
#catch e
#    println(e)
#    Pkg.add("FastaIO")
#    using FastaIO
#end

# using ArgParse
#try
#    Pkg.installed("ArgParse")
#    using ArgParse
#catch e
#    println(e)
#    Pkg.add("ArgParse")
#    using ArgParse
#end


function parse_commandline()

    settings = ArgParseSettings("Count index frequencies in a fastq file",
                                commands_are_required = false,
                                version = "1.0",
                                add_version = true)

    @add_arg_table settings begin
        "input_fasta_file"
            help = "input unfiltred fasta file"
            required = true
        "--prefix", "-p"
            help = "output prefix"
            arg_type = String
        "--length_cutoff", "-l"
            help = "input minimum length_cutoff to filter the raw fasta file"
            arg_type = Int64
            default = 1000
        "--reverse", "-r"
            help = "reverse the selection, write sequences shorter than this cutoff?"
            action = :store_true
        "--force", "-f"
            help = "force to remove existing output file?"
            action = :store_true
    end

    return parse_args(settings)
end


function filter_fasta_by_length(input_fasta_file, length_cutoff, output_file, reverse=false)
    # expand filepath
    output_file = expanduser(output_file)
    input_fasta_file = expanduser(input_fasta_file)

    # output handle
    FastaWriter(output_file) do output_handle
        # fasta iterator 
        FastaReader(input_fasta_file) do fasta_records
            for (desc, seq) in fasta_records
                #println("$desc : $seq")
                if reverse
                    if length(seq) < length_cutoff
                        writeentry(output_handle, desc, seq)
                    end
                else
                    if length(seq) >= length_cutoff
                    writeentry(output_handle, desc, seq) 
                    end
                end
            end
        end
    end  
end



# -----------
#  main
# -----------

function main()
    args = parse_commandline()
    #println("Parsed args:")
    #for (arg,val) in args
    #    println("  $arg  =>  $val")
    #end

    # output_file
    length_cutoff = args["length_cutoff"]
    base_name = basename(args["input_fasta_file"])
    file_stem = splitext(base_name)[1]
    if args["prefix"] != nothing 
        prefix = args["prefix"]
    else
        prefix = file_stem * "_min" * "$length_cutoff"
    end
    output_file = string(prefix, ".fa")

    # In case of file exists
    if isfile(output_file)
        if ! args["force"]
           @error "$output_file exists, please rename it and re-run again."
           exit()
        else
            @warn "$output_file exists, will be overwritten."
        end   
    end 

    # filter fasta by length_cutoff
    filter_fasta_by_length(args["input_fasta_file"], length_cutoff, output_file, args["reverse"])

end


# call main
main()

