function calculateAtomsPerReaction(reactions,data_dictionary, local_atom_dictionary)
    list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
    number_of_metabolites = length(list_of_metabolite_symbols_model)
    tmp_array::Array{AbstractString} = AbstractString[]
    atom_array = zeros(number_of_metabolites,6)
    for line in reactions
        line=replace(line, r":+", ",")
        #@show line
        split_line = split(line, ",")
        if(size(split_line,1)>2)
            rxn_id=split_line[1]
            name = split_line[2]
            flux = split_line[4]
            rxn = split_line[3]
        else
            name = split_line[1]
            rxn = split_line[2]
        end
        split_rxn = split(strip(rxn), "-->")
        #@show split_rxn
        total_atoms = zeros(6)
        for idx in (1,size(split_rxn,1))
            section = split_rxn[idx]
            #@show section
            species = split(strip(section), "+")
            for item in species
                #@show item
                coeff = 1
                #if we have a coefficent other than 1
                if(searchindex(item, '*')!=0)
                    staridx = searchindex(item, '*')
                    coeff = parse(item[1:staridx-1])
                    item = item[staridx+1:end]
                end
                atoms = local_atom_dictionary[item]
                if(idx == 1) #reactants
                    total_atoms = atoms*coeff+total_atoms
                else #products
                    total_atoms = -atoms*coeff+total_atoms
                end
            end
        end
        if(maximum(abs(total_atoms))>0) #show unbalanced reactions
            @show name
            @show total_atoms
        end

    end
end
function generate_local_atom_arr(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})
  # how many metabolite symbols do we have in *the model*?
  list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
  number_of_metabolites = length(list_of_metabolite_symbols_model)
  # initialize -
  tmp_array::Array{AbstractString} = AbstractString[]
  atom_array = zeros(number_of_metabolites,6)
  local_dictionary::Dict{AbstractString,Any} = Dict{AbstractString,Any}()
  # load the atom file -
  try
    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)
          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end
    # ok, create a local dictionary w/the atom records -
    for record in tmp_array
      # split -
      split_array = split(record,",")
      # get my key -
      key = split_array[1]  # Metabolite symbol -
      # local array -
      local_atom_array = zeros(6)
      local_atom_array[1] = parse(Float64,split_array[2]) # C
      local_atom_array[2] = parse(Float64,split_array[3]) # H
      local_atom_array[3] = parse(Float64,split_array[4]) # N
      local_atom_array[4] = parse(Float64,split_array[5]) # O
      local_atom_array[5] = parse(Float64,split_array[6]) # P
      local_atom_array[6] = parse(Float64,split_array[7]) # S
      # store -
      local_dictionary[key] = local_atom_array
    end
      catch err
    showerror(STDOUT, err, backtrace());println()
  end
    return local_dictionary
end
function checkAllBalances()
    path_to_atom_file = "Atom.txt"
    data_dictionary =DataDictionary(0,0,0)
    local_atom_array=generate_local_atom_arr(path_to_atom_file, data_dictionary)
    calculateAtomsPerReaction(data_dictionary["list_of_reaction_strings"],data_dictionary, local_atom_array)
end
