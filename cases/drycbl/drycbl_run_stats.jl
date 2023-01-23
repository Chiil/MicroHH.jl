## Loading packages.
using Statistics
# using Profile
using MicroHH
using MicroHH.Stats


## Loading settings and stats function.
include("drycbl_settings.jl")


## Initialize the model.
n_domains = 1
m = Model("drycbl", n_domains, settings)


## Load the restart data.
load_model!(m)


## Prepare stats file and create stats function with dimensions needed.
# let
#     fid = open_stats("drycbl_stats.h5")
#     g = m.grid[1]
#     add_dimension!(fid, "x", g.x[g.is:g.ie])
#     add_dimension!(fid, "y", g.y[g.js:g.je])
#     add_dimension!(fid, "z", g.z[g.ks:g.ke])
#     add_dimension!(fid, "zh", g.zh[g.ks:g.keh])
#     close(fid)
# end

let
    fid = open_stats("drycbl_cross.h5")
    g = m.grid[1]
    add_dimension!(fid, "x", g.x[g.is:g.ie])
    add_dimension!(fid, "xh", g.xh[g.is:g.ie+1])
    add_dimension!(fid, "z", g.z[g.ks:g.ke])
    add_dimension!(fid, "zh", g.zh[g.ks:g.keh])
    close(fid)
end

function do_stats(m::Model)
    g = m.grid[1]; f = m.fields[1]; t = m.timeloop[1]
    fid = open_stats("drycbl_stats.h5")

    # 1. Add the time
    add_time_record!(fid, m.timeloop[1].time, t.iter)

    # 2. Add desired cross sctions.
    s = @view f.s[g.is:g.ie, g.js, g.ks:g.ke]
    add_record!(fid, s, "s_xz", ("x", "z"), t.iter)

    # 3. Add profile
    u = @view f.u[g.is:g.ie, g.js:g.je, g.ks:g.ke ]
    v = @view f.v[g.is:g.ie, g.js:g.je, g.ks:g.ke ]
    w = @view f.w[g.is:g.ie, g.js:g.je, g.ks:g.keh]
    s = @view f.s[g.is:g.ie, g.js:g.je, g.ks:g.ke ]

    add_record!(fid, mean(u, dims=(1, 2))[:], "u", ("z", ), t.iter)
    add_record!(fid, mean(v, dims=(1, 2))[:], "v", ("z", ), t.iter)
    add_record!(fid, mean(w, dims=(1, 2))[:], "w", ("zh",), t.iter)
    add_record!(fid, mean(s, dims=(1, 2))[:], "s", ("z", ), t.iter)

    add_record!(fid, var(u, dims=(1, 2))[:], "u_2", ("z", ), t.iter)
    add_record!(fid, var(v, dims=(1, 2))[:], "v_2", ("z", ), t.iter)
    add_record!(fid, var(w, dims=(1, 2))[:], "w_2", ("zh",), t.iter)
    add_record!(fid, var(s, dims=(1, 2))[:], "s_2", ("z", ), t.iter)

    close(fid)
end

function do_cross(m::Model)
    g = m.grid[1]; f = m.fields[1]; t = m.timeloop[1]
    fid = open_stats("drycbl_cross.h5")

    # 1. Add the time
    add_time_record!(fid, m.timeloop[1].time, t.iter)

    # 2. Add desired cross sctions.
    s = @view f.s[g.is:g.ie, g.js, g.ks:g.ke]
    add_record!(fid, s, "s_xz", ("x", "z"), t.iter)

    u = @view f.u[g.is:g.ie+1, g.js, g.ks:g.ke]
    add_record!(fid, u, "u_xz", ("xh", "z"), t.iter)

    w = @view f.w[g.is:g.ie, g.js, g.ks:g.keh]
    add_record!(fid, w, "w_xz", ("x", "zh"), t.iter)

    p = @view f.p[g.is:g.ie, g.js, g.ks:g.ke]
    add_record!(fid, p, "p_xz", ("x", "z"), t.iter)

    close(fid)
end


## Run the model.
in_progress = prepare_model!(m)
# do_stats(m)
do_cross(m)

while in_progress
    global in_progress = step_model!(m)
    # @profile global in_progress = step_model!(m)

    # if mod(round(Int, m.timeloop[1].time), 60) == 0
    #     do_stats(m)
    # end

    if mod(round(Int, m.timeloop[1].time), 10) == 0
        do_cross(m)
    end
end


## Closing actions.
output_timer_model!(m)
