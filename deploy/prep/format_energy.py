
hf = open("tmp_data/E_HF.dat", 'r')
ci = open("tmp_data/E_CI.dat", 'r')

hf_line = hf.readline()
ci_line = ci.readline()

hf.close()
ci.close()

hf_new = open("E_HF.dat", 'w')
ci_new = open("E_CI.dat", 'w')

hf_new.write(hf_line.split()[4] + "\n")
ci_new.write(ci_line.split()[3] + "\n")

hf_new.close()
ci_new.close()