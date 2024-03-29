using .StencilBuilder


const ifactor = 1_000_000


mutable struct Timeloop
    # Float numbers.
    start_time::Float64
    end_time::Float64
    save_time::Float64
    check_time::Float64
    dt::Float64
    time::Float64

    # Integer values.
    istart_time::Int64
    iend_time::Int64
    isave_time::Int64
    icheck_time::Int64
    idt::Int64
    itime::Int64

    rkstep::Int64
    iter::Int64
end


function Timeloop(settings::Dict)
    d = settings["timeloop"]

    start_time = d["start_time"]
    end_time = d["end_time"]
    save_time = d["save_time"]
    check_time = d["check_time"]
    dt = d["dt"]
    time = start_time

    istart_time = round(Int64, start_time * ifactor)
    iend_time = round(Int64, end_time * ifactor)
    isave_time = round(Int64, save_time * ifactor)
    icheck_time = round(Int64, check_time * ifactor)
    idt = round(Int64, dt * ifactor)
    itime = istart_time

    rkstep = 1
    iter = istart_time ÷ idt

    Timeloop(
        start_time, end_time, save_time, check_time, dt, time,
        istart_time, iend_time, isave_time, icheck_time, idt, itime,
        rkstep, iter)
end


function integrate_time_kernel!(
    a, at,
    rkstep, dt,
    is, ie, js, je, ks, ke)

    c_a = [0, -5//9, -153//128];
    c_b = [1//3, 15//16, 8//15];

    rkstep_next = (rkstep + 1) == 4 ? 1 : rkstep + 1

    c_b_step = c_b[rkstep] * dt
    c_a_step = c_a[rkstep_next]

    @fast3d begin
        @fd (a[i, j, k], at[i, j, k]) a += c_b_step*at
        @fd (a[i, j, k], at[i, j, k]) at *= c_a_step
    end
end


function integrate_time!(
    f::Fields, g::Grid, t::Timeloop)

    # Make sure the dt is in Float32 if the array is.
    dt = convert(eltype(f.u), t.dt)

    integrate_time_kernel!(
        f.u, f.u_tend,
        t.rkstep, dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    integrate_time_kernel!(
        f.v, f.v_tend,
        t.rkstep, dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    integrate_time_kernel!(
        f.w, f.w_tend,
        t.rkstep, dt,
        g.is, g.ie, g.js, g.je, g.ks, g.keh)

    integrate_time_kernel!(
        f.s, f.s_tend,
        t.rkstep, dt,
        g.is, g.ie, g.js, g.je, g.ks, g.ke)

    for scalar_name in keys(f.scalars)
        integrate_time_kernel!(
            f.scalars[scalar_name], f.scalars_tend[scalar_name],
            t.rkstep, dt,
            g.is, g.ie, g.js, g.je, g.ks, g.ke)
    end
end


function step_time!(t::Timeloop)
    t.rkstep += 1

    if t.rkstep > 3
        t.rkstep = 1
        t.itime += t.idt
        t.time = convert(Float64, t.itime) / ifactor
        t.iter += 1
    end
end


function get_sub_dt(t::Timeloop)
    c_b = [1//3, 15//16, 8//15];
    return c_b[t.rkstep]*t.dt;
end


function is_in_progress(t::Timeloop)
    return t.itime < t.iend_time
end


function is_save_time(t::Timeloop)
    return t.itime % t.isave_time == 0 && t.rkstep == 1 && t.itime != t.istart_time
end


function is_check_time(t::Timeloop)
    return t.itime % t.icheck_time == 0
end
