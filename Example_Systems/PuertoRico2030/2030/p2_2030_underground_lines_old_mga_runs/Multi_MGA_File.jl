"""
This script compiles capacity and network results for parallel MGA runs into two files for ease of processing.

*** Modified Version for multiple parallel MGA runs ***

"""

using CSV
using DataFrames
using YAML


# Folder Path - Input path to caserunner here
CaseRunner_Path = raw"/home/aerivera/SimpleGenXCaseRunner_Local"

yml2dict(path::AbstractString) = YAML.load_file(path)
dict2yml(d::Dict, path::AbstractString) = YAML.write_file(path, d)

csv2dataframe(path::AbstractString) = CSV.read(path, header=1, DataFrame)
dataframe2csv(df::DataFrame, path::AbstractString) = CSV.write(path, df)

"""
Function to add case label column
"""

function label_col(df::DataFrame, MGA_Run::AbstractString)
    row_num = nrow(df)
    A = Vector{AbstractString}(undef, row_num)
    for i in 1:row_num
        A[i] = MGA_Run
    end
    df.MGA_Iteration = A
    return df
end

function make_label(case_num::Int64, direction::AbstractString, run_num::Int64)
    case = "Case_"
    c_num = string(case_num)
    r_num = string(run_num)
    label = case*c_num*"_"*direction*"_"*r_num
    return label
end

function create_csvs(net_df::DataFrame, cap_df::DataFrame, path::AbstractString)
    out_capacity = "Raw_Capacity.csv"
    out_network = "Raw_Network.csv"

    cap_path = joinpath(path, out_capacity)
    net_path = joinpath(path, out_network)
    
    dataframe2csv(cap_df, cap_path)
    dataframe2csv(net_df, net_path)
end


function compile_dfs(all_cases_path::AbstractString)

    opt = "Optimal Solution"
    max = "MGAResults_max"
    up = "max"
    min = "MGAResults_min"
    down = "min"
    results = "Results"
    capacity = "capacity.csv"
    network = "network_expansion.csv"

    all_entries = readdir(all_cases_path, join=true)
    counter = 0
    for f in all_entries # entering Cases folder
        counter += 1
        indiv_case_files = readdir(f, join=false) # Reading all cases inside (Case_1, 2 etc)
        if counter == 1
            for k in indiv_case_files # Separate for loop here to ensure construction of dataframes
                if k == results
                    k_path = joinpath(f, k)
                    opt_csvs = readdir(k_path, join=false)
                    for m in opt_csvs # Getting optimal solution network and capacity first
                        if m == capacity
                            df_path = joinpath(k_path, m)
                            raw_cap_df = csv2dataframe(df_path)
                            global raw_cap_df = label_col(raw_cap_df, opt)
                        elseif m == network
                            df_path = joinpath(k_path, m)
                            raw_network_df = csv2dataframe(df_path)
                            global raw_network_df = label_col(raw_network_df, opt)
                        end
                    end
                end   
            end
        end
        for k in indiv_case_files
            if k == max || k == min
                k_path = joinpath(f, k)
                runs = readdir(k_path, join=false)
                run_num = 0
                for l in runs
                    run_num += 1
                    if k == max
                        l_path = joinpath(k_path, l)
                        max_csvs = readdir(l_path, join=false)
                        for m in max_csvs # Getting optimal solution network and capacity first
                            if m == capacity
                                df_path = joinpath(l_path, m)
                                new_cap_df = csv2dataframe(df_path)
                                label = make_label(counter, up, run_num)
                                new_cap_df = label_col(new_cap_df, label)
                                append!(raw_cap_df, new_cap_df)
                                
                            elseif m == network
                                df_path = joinpath(l_path, m)
                                new_network_df = csv2dataframe(df_path)
                                label = make_label(counter, up, run_num)
                                new_network_df = label_col(new_network_df, label)
                                append!(raw_network_df, new_network_df)
                            end
                        end
                    elseif k == min
                        l_path = joinpath(k_path, l)
                        min_csvs = readdir(l_path, join=false)
                        for m in min_csvs # Getting optimal solution network and capacity first
                            if m == capacity
                                df_path = joinpath(l_path, m)
                                new_cap_df = csv2dataframe(df_path)
                                label = make_label(counter, down, run_num)
                                new_cap_df = label_col(new_cap_df, label)
                                append!(raw_cap_df, new_cap_df)
                            elseif m == network
                                df_path = joinpath(l_path, m)
                                new_network_df = csv2dataframe(df_path)
                                label = make_label(counter, down, run_num)
                                new_network_df = label_col(new_network_df, label)
                                append!(raw_network_df, new_network_df)
                            end
                        end
                    end
                end
            end
        end
    end
    return raw_network_df, raw_cap_df
end

# Moving all results files into one consolidated file
function main()
    present_dir = pwd()
    Cases = "Cases"
    Comp_Results = "MGA_Results"
    mkdir(Comp_Results)
    Final_path = joinpath(present_dir, Comp_Results)
    all_cases_path = joinpath(CaseRunner_Path, Cases)

    raw_network_df, raw_cap_df = compile_dfs(all_cases_path)
    create_csvs(raw_network_df, raw_cap_df, Final_path)
end 

main()