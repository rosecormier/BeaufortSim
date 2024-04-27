function ReadInputFile(filename)
    
    input_params = Dict()
    
    for (i, line) in enumerate(eachline(filename))
        
        var_len, line_len = 0, length(line)
        if (line_len > 1 && line[1] != '#')

            for char in line
                if char == '='
                    break
                else
                    var_len += 1
                end
            end

            var = strip(line[1:var_len])
            line = strip(line[var_len:end])

            val_len = 0
            for char in line
                if char == '#' 
                    break
                else
                    val_len += 1
                end
            end

            val_str = strip(line[2:val_len]) #Index 1 is "="

            val = tryparse(Int, val_str)
            if isnothing(val)
                val = parse(Float64, val_str)
            end

            input_params[var] = val

        end    
    end
    
    input_params
    
end