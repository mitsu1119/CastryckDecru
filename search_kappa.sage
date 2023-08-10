import sys

load("attack.sage")

file_name = sys.argv[1]

kappa_choice_params = load(file_name + ".sobj")
res = choiced_kappa(kappa_choice_params)

save(res, file_name)
