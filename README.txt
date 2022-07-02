Autonomous hybrid green power system consists of PV, WT, BT, SC, HT, EL, and FC.
P_PV.m, P_WT.m, E_BT_fun.m, E_SC_fun.m, E_H2_fun.m, V_EL.m, and V_FC.m
models the equipments respectively.
start the program by reading Main.m to get the gist of it.
MegaRun_over_Main.m is used to run Main.m for numerous times
(for testing or simulation purposes).
details.m simulates the energy management system (EMS).
details.m (EMS) is based on three fuzzy logic controllers (FLC) proposed by R. Zahedi.
R. Zahedi: https://doi.org/10.1016/j.energy.2020.117935.
OP stands for Optimal Parameters in fuzzy logic controller.
FLC_optimized produces the control signals for the storage system based on
the parameters of the fuzzy logic
membership function and the rules given to it.
OR stands for Optimal Rules in fuzzy logic controller.
FixPosition.m makes sure that the parameters of the membership function
of FLC are within their limits
through out the optimization process.
OS stands for Optimal Size in the optimal design process.
FixSize.m makes sure that the size of the equipments stay within their limits
through out the optimization process.
CostFunction.m caculated the cost (initial cost + replacement cost +
operation and maintenance + degradation penalty)
for parPSO and parEPSO programs.
Through Constraints.m, the functional limitations of the renewable energy system
is passed to the optimization algorithms.
parPSO and parEPSO are PSO and E-PSO algorithms implemented using parfor-loop
as in parallel computation.
To limit the calculations, the hierarchical optimization is applied.

PV_test.xlsx and WT_test.xlsx contain test data for PV and WT for
a 24 hr (24 data point) beta test.
test_data.xlsx contains the load, irradiation, ambient temperature,
and wind speed for a 24 hr beta test.
run_time.xlsx records the time that the simulation ends.
datetime_value.m returns the date and time of the moment it gets called
as a 14 digit integer number.
details_FLC_log.xlsx records the optimal results found at the end of
each optimization process for FLC and OS.
