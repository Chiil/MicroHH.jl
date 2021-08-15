using .StencilBuilder

mutable struct Timeloop
    start_time::Float64
    end_time::Float64
    dt::Float64

    time::Float64
    rkstep::Int
end

function Timeloop(d::Dict)
    start_time = d["start_time"]
    end_time = d["end_time"]
    dt = d["dt"]
    rkstep = 1

    Timeloop(start_time, end_time, dt, start_time, rkstep)
end

function integrate_time_kernel!(
    a, at,
    rkstep, dt,
    is, ie, js, je, ks, ke)

    c_a = [0, 5/9, -153/128];
    c_b = [1/3, 15/16, 8/15];

    rkstep_next = (rkstep + 1) == 4 ? 1 : rkstep + 1

    c_b_step = c_b[rkstep] * dt
    c_a_step = c_a[rkstep_next]

    @fast3d begin
        @fd (a, at) a += c_b_step*at
        @fd (a, at) at *= c_a_step
    end
end

function integrate_time!(
    fields::Fields, grid::Grid, timeloop::Timeloop)

    integrate_time_kernel!(
        fields.u, fields.u_tend,
        timeloop.rkstep, timeloop.dt,
        grid.is, grid.ie, grid.js, grid.je, grid.ks, grid.ke)

    integrate_time_kernel!(
        fields.v, fields.v_tend,
        timeloop.rkstep, timeloop.dt,
        grid.is, grid.ie, grid.js, grid.je, grid.ks, grid.ke)

    integrate_time_kernel!(
        fields.w, fields.w_tend,
        timeloop.rkstep, timeloop.dt,
        grid.is, grid.ie, grid.js, grid.je, grid.ks, grid.keh)

    integrate_time_kernel!(
        fields.s, fields.s_tend,
        timeloop.rkstep, timeloop.dt,
        grid.is, grid.ie, grid.js, grid.je, grid.ks, grid.ke)
end

function step_time!(timeloop::Timeloop)
    timeloop.rkstep += 1

    if timeloop.rkstep > 3
        timeloop.rkstep = 1
        timeloop.time += timeloop.dt
    end
end

function get_sub_dt(t::Timeloop)
    c_b = [1/3, 15/16, 8/15];
    return c_b[t.rkstep]*t.dt;
end
