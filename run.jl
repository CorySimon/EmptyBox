using DataFrames

# Simulate background field model with different background energies
U0 = linspace(-20*1000/8.314+298,0, 50)
n = length(U0)
for i = 1:n
    @printf("Simulating %d/%d\n", i, n)
    E = U0[i]
    run(`./simulate $E 20000`)
end

# collect results in a data frame
df = DataFrame(U0=U0, L1=zeros(n), L5pt8=zeros(n), L35=zeros(n), L65=zeros(n))
for i = 1:length(U0)
    df_temp = readtable(@sprintf("results/SIM_U%.2f.txt", U0[i]), skipstart=2)
    df[:L1][i] = df_temp[:Loading_vSTP_v_][1]
    df[:L5pt8][i] = df_temp[:Loading_vSTP_v_][2]
    df[:L35][i] = df_temp[:Loading_vSTP_v_][3]
    df[:L65][i] = df_temp[:Loading_vSTP_v_][4]
end
writetable("loading_of_U.csv", df)
