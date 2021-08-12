mutable struct Timeloop
    start_time::Float64
    end_time::Float64
    dt::Float64

    time::Float64
end

function Timeloop(d::Dict)
    start_time = d["start_time"]
    end_time = d["end_time"]
    dt = d["dt"]

    Timeloop(start_time, end_time, dt, start_time)
end

function step_time!(timeloop::Timeloop)
    timeloop.time += timeloop.dt
    return timeloop.time < timeloop.end_time
end
