import spectrum

filepath=r"C:\\Users\\ALber\\Desktop\\Inp2.spe"
list_erg,list_eff=spectrum.read(filepath)

print(str(spectrum.areaTPA(list_erg,list_eff,peak=305,avg_number=5,width_base=5,base_method=0)))  