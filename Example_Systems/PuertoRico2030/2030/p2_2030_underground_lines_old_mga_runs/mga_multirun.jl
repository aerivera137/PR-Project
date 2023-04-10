"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

"""
This script allows MGA runs to be parallelized using system at https://github.com/cfe316/SimpleGenXCaseRunner.git. Please ensure the joblocation in caserunner.jl is set to "BATCH"
To use it, Change CaseRunner_Path to your local caserunner folder. Move this file into your GenX Running Folder, like "..GenX/Example_Systems/RealSystemExample/MGA_ISONE_Trizone_FullTimeseries/".
Set everything up as though this were a normal MGA run, then specify number of MGA sessions you want to run in parallel and run this file.
"""


using CSV
using DataFrames
using YAML


# Folder Path - Input path to caserunner here
CaseRunner_Path = raw"/home/aerivera/SimpleGenXCaseRunner_Local"
Myopic_Path = raw"/home/aerivera/GenX-Local/GenX/Example_Systems/PuertoRico2030"

yml2dict(path::AbstractString) = YAML.load_file(path)
dict2yml(d::Dict, path::AbstractString) = YAML.write_file(path, d)

csv2dataframe(path::AbstractString) = CSV.read(path, header=1, DataFrame)
dataframe2csv(df::DataFrame, path::AbstractString) = CSV.write(path, df)

"""
Functions to modifying CSVs in pwd
"""

# Replaces one value in a dataframe with a special key

function insert_special_key_csv(insert_loc::DataFrame, replacements::DataFrame, key::AbstractString)
    println("in insert special key")
    replacement = "__SPECIAL_"*key*"__"
    element = insert_loc[4,4] # Change at will to something that works, this is the value of the cell we're replacing
    println(element)
    new_col(replacements, element, key)
    insert_loc[!, 4] = convert(Vector{Union{Float64, String}}, insert_loc[:, 4])
    insert_loc[4, 4] = replacement
end

"""
For a given change number, value, and dataframe, form necessary new rows and columns for replacements.csv
"""

function add_case(df::DataFrame)
    r_num = nrow(df)
    println(r_num)
    println("Adding Case")
    push!(df, df[r_num,:])
    df[r_num + 1, 1] = r_num + 1
    df[r_num + 1, 2] = string(r_num)
    return df
end

function new_col(df::DataFrame, value::Float64, key::AbstractString)
    Case = "Case"
    c_new = DataFrame(Case => 1, key => value)
    leftjoin!(df, c_new, on =:Case)
    return df
end

function initialize_df()
    df = DataFrame(Case = 1, Notes = "First")
    return df
end

function initialize_repcsv(df::DataFrame, rep_path::AbstractString)
    replacement = "replacements.csv"
    path = joinpath(rep_path, replacement)
    println(path)
    println(df)
    CSV.write(path, df)
end

# Copying all files for MGA run into temp_dir, the path for template folder
function copy_files(temp_dir::AbstractString, myopic_path::AbstractString)
    cp(myopic_path, temp_dir, force = true)
end


# Main executable
function main()
    # Path strings
    template = "template"
    settings = "Settings"
    genx_settings = "genx_settings.yml"
    execute = "caserunner.jl"
    key = "key"
    inpath = pwd()
    #infile = "Inputs"
    #infile = "p2_2030_underground_lines"
    #randfile_path = joinpath(inpath, infile,"Fuels_data.csv")
    randfile_path = joinpath(inpath, "Fuels_data.csv")
    setup_path = joinpath(inpath, settings, genx_settings)
    temp_path = joinpath(CaseRunner_Path, template)
    println("done reading paths")

    s = yml2dict(setup_path)
    replace_df = initialize_df()
    num_iters = s["ModelingToGenerateAlternativesParallel"]
    num_iters = convert(Int64, num_iters)
    read_file = csv2dataframe(randfile_path)
    println(randfile_path)
    insert_special_key_csv(read_file, replace_df, key)
    for i in 2:num_iters
        add_case(replace_df)
    end
    dataframe2csv(read_file, randfile_path)


    # Move things to where they ought to be and instantiate 
    copy_files(temp_path, Myopic_Path)
    initialize_repcsv(replace_df, CaseRunner_Path)
    

    # Run caserunner
    cd(CaseRunner_Path)
    execute_path = joinpath(CaseRunner_Path, execute)
    include(execute_path)

    # TODO - Include Post Process
end


main()