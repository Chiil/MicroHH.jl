## The cells in this file can only be run after running those in drycbl_run_plot.jl


## Transfer the state from the high to the low resolution and plot it.
transfer_state(f1, f2, g1, g2)
node1[] = s_bot1 .- mean(s_bot1)
node2[] = s_bot2 .- mean(s_bot2)


## Run the model for another 900 sec.
transfer_state(f1, f2, g1, g2)
m.timeloop[1].end_time += 900.
m.timeloop[1].iend_time += 900 * MicroHH.ifactor
in_progress = true
while in_progress
    global in_progress = step_model!(m)
    node1[] = s_bot1 .- mean(s_bot1)
    node2[] = s_bot2 .- mean(s_bot2)
end
