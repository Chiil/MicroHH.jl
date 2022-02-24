## Loading packages.
using Statistics
using MicroHH
# using Profile


## Loading settings and stats function.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)


## Load the restart data.
load_model!(m)


## Prepare stats file and create stats function.
let
    fid = open_stats("drycbl_stats.h5")
    g = m.grid[1]
    add_dimension!(fid, "x", g.x[g.is:g.ie])
    add_dimension!(fid, "y", g.y[g.js:g.je])
    close(fid)
end

function do_stats(m::Model)
    g = m.grid[1]; f = m.fields[1]; t = m.timeloop[1]
    fid = open_stats("drycbl_stats.h5")
    add_time_record!(fid, m.timeloop[1].time, t.iter)
    s = @view f.s[g.is:g.ie, g.js:g.je, g.ks]
    add_record!(fid, s .- mean(s), "s_xy", ("x", "y"), t.iter)
    close(fid)
end


## Run the model.
in_progress = prepare_model!(m)
do_stats(m)

while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)

    do_stats(m)
end


## Closing actions.
output_timer_model!(m)
