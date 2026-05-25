function _nonempty_tokens(line)
    return split(strip(line))
end

function _parse_openfast_output_row(line, ncolumns)
    tokens = _nonempty_tokens(line)
    length(tokens) == ncolumns || return nothing
    values = Vector{Float64}(undef, ncolumns)
    for (idx, token) in enumerate(tokens)
        value = tryparse(Float64, token)
        isnothing(value) && return nothing
        values[idx] = value
    end
    return values
end

"""
    readOpenFASTOutputTable(filename)

Read a text OpenFAST, AeroDyn-driver, or module-driver tabular output file.
The returned named tuple contains `channels`, `units`, and a dense numeric
`data` matrix. Non-tabular preamble lines are skipped, embedded NUL bytes are
treated as spaces, and malformed trailing lines are ignored.
"""
function readOpenFASTOutputTable(filename)
    isfile(filename) || throw(ArgumentError("OpenFAST output file does not exist: $filename"))

    text = replace(read(filename, String), Char(0) => ' ')
    lines = split(text, '\n')
    header_index = findfirst(line -> startswith(strip(line), "Time"), lines)
    isnothing(header_index) &&
        throw(ArgumentError("OpenFAST output table header starting with Time was not found"))

    channels = _nonempty_tokens(lines[header_index])
    isempty(channels) && throw(ArgumentError("OpenFAST output table has no channel names"))
    units = header_index < length(lines) ? _nonempty_tokens(lines[header_index + 1]) : String[]

    rows = Vector{Vector{Float64}}()
    for line in lines[(header_index + 2):end]
        values = _parse_openfast_output_row(line, length(channels))
        isnothing(values) || push!(rows, values)
    end
    isempty(rows) && throw(ArgumentError("OpenFAST output table has no numeric rows"))

    data = Matrix{Float64}(undef, length(rows), length(channels))
    for (row_index, values) in enumerate(rows)
        data[row_index, :] .= values
    end

    return (channels = channels, units = units, data = data)
end

function _openfast_table(table_or_filename)
    return table_or_filename isa AbstractString ?
           readOpenFASTOutputTable(table_or_filename) :
           table_or_filename
end

function _openfast_channel_index(table, channel)
    index = findfirst(==(channel), table.channels)
    isnothing(index) && throw(ArgumentError("channel $(repr(channel)) was not found"))
    return index
end

function _validated_openfast_metric_rows(reference, candidate, rows)
    max_rows = min(size(reference.data, 1), size(candidate.data, 1))
    if isnothing(rows)
        return 1:max_rows
    end
    selected_rows = collect(rows)
    isempty(selected_rows) && throw(ArgumentError("rows must not be empty"))
    all(row -> row isa Integer && 1 <= row <= max_rows, selected_rows) ||
        throw(ArgumentError("rows must be integer indices in the common output-table range"))
    return selected_rows
end

"""
    openfastOutputChannelMetrics(reference, candidate, channels; rows=nothing)

Compare selected output `channels` between two OpenFAST-style output tables or
filenames. `candidate - reference` errors are evaluated over `rows`, or over
the common row range when `rows` is omitted. The return value is a dictionary
from channel name to metrics: `n`, `rmse`, `mean_bias`, `mean_abs_error`,
`max_abs_error`, final values, and final error.
"""
function openfastOutputChannelMetrics(reference, candidate, channels; rows = nothing)
    reference_table = _openfast_table(reference)
    candidate_table = _openfast_table(candidate)
    selected_rows = _validated_openfast_metric_rows(reference_table, candidate_table, rows)

    metrics = Dict{String,NamedTuple}()
    for channel in channels
        channel_name = String(channel)
        reference_index = _openfast_channel_index(reference_table, channel_name)
        candidate_index = _openfast_channel_index(candidate_table, channel_name)
        reference_values = reference_table.data[selected_rows, reference_index]
        candidate_values = candidate_table.data[selected_rows, candidate_index]
        error = candidate_values .- reference_values
        abs_error = abs.(error)
        n = length(error)
        metrics[channel_name] = (
            n = n,
            rmse = sqrt(sum(error .^ 2) / n),
            mean_bias = sum(error) / n,
            mean_abs_error = sum(abs_error) / n,
            max_abs_error = maximum(abs_error),
            reference_final = reference_values[end],
            candidate_final = candidate_values[end],
            final_error = error[end],
        )
    end
    return metrics
end
