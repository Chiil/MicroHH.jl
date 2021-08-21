using .StencilBuilder

mutable struct Timeloop
    start_time::Float64
    end_time::Float64
    save_time::Float64
    check_time::Float64
    dt::Float64

    time::Float64
    rkstep::Int
end

function Timeloop(d::Dict)
    start_time = d["start_time"]
    end_time = d["end_time"]
    save_time = d["save_time"]
    check_time = d["check_time"]
    dt = d["dt"]
    rkstep = 1
    time = start_time

    Timeloop(start_time, end_time, save_time, check_time, dt, time, rkstep)
end

function integrate_time_kernel!(
    a, at,
    rkstep, dt,
    is, ie, js, je, ks, ke)

    c_a = [0, 5//9, -153//128];
    c_b = [1//3, 15//16, 8//15];

    rkstep_next = (rkstep + 1) == 4 ? 1 : rkstep + 1

    c_b_step = c_b[rkstep] * dt
    c_a_step = c_a[rkstep_next]

    @fast3d begin
        @fd (a, at) a += c_b_step*at
        @fd (a, at) at *= c_a_step
    end
end

function integrate_time!(
    f::Fields, g::Grid, t::Timeloop)

    integrate_time_kernel!(
        f.u, f.u_tend,
        t.rkstep, t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    integrate_time_kernel!(
        f.v, f.v_tend,
        t.rkstep, t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    integrate_time_kernel!(
        f.w, f.w_tend,
        t.rkstep, t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.keh)

    integrate_time_kernel!(
        f.s, f.s_tend,
        t.rkstep, t.dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)
end

function step_time!(t::Timeloop)
    t.rkstep += 1

    if t.rkstep > 3
        t.rkstep = 1
        t.time += t.dt
    end
end

function get_sub_dt(t::Timeloop)
    c_b = [1/3, 15/16, 8/15];
    return c_b[t.rkstep]*t.dt;
end
