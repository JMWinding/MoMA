gold:
cdmalen = goldlen
PD = prelen * cdmalen * T (each chip of gold repeated prelen times)
SI = cdmalen * T

goldman: 
cdmalen = 2 * goldlen
PD = prelen * cdmalen * T (each chip of goldman repeated prelen times)
SI = cdmalen * T

plain: first half chips are all "1"s for bit "1", and no pulse for bit "-1"
cdmalen = 2 * goldlen
PD = prelen * cdmalen * T (each chip of goldman repeated prelen times)
SI = cdmalen * T

man: first half chips are all "1"s for bit "1", and second half chips are all "1"s for bit "-1"
cdmalen = 2 * goldlen
PD = prelen * cdmalen * T (each chip of goldman repeated prelen times)
SI = cdmalen * T

plain0: long symbol interval with only one pulse in the beginning, and no pulse for bit "-1"
cdmalen = 2 * goldlen
PD = prelen * cdmalen * T (goldman repeated prelen times)
SI = 2 * T