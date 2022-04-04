## Loading packages.
using MicroHH
using MicroHH.Stats
using Statistics
# using Profile


## Loading settings.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings, float_type)


## Load the restart data.
load_model!(m)


## Create the cross section file and function.
let
    fid = open_stats("drycbl_cross.h5")
    g = m.grid[1]
    add_dimension!(fid, "x", g.x[g.is:g.ie])
    add_dimension!(fid, "z", g.z[g.ks:g.ke])
    add_dimension!(fid, "xh", g.xh[g.is:g.ie])
    add_dimension!(fid, "zh", g.zh[g.ks:g.keh])
    close(fid)
end

function do_cross(m::Model)
    g = m.grid[1]; f = m.fields[1]; t = m.timeloop[1]
    fid = open_stats("drycbl_cross.h5")

    # 1. Add the time
    add_time_record!(fid, m.timeloop[1].time, t.iter)

    # 2. Add desired cross sections.
    s = @view f.s[g.is:g.ie, g.js, g.ks:g.ke]
    add_record!(fid, s, "s_xz", ("x", "z"), t.iter)

    close(fid)
end


## Set up the model.
in_progress = prepare_model!(m)
do_cross(m)


## Run the model.
while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)

    if mod(m.timeloop[1].itime, 20_000) == 0
        do_cross(m)
    end
end
