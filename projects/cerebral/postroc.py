
delta=0.01

def beolvaso(filenev): 

    fbe=open(filenev,"r")
    l=fbe.readlines()
    for i in range(len(l)):
        l[i]=l[i][:-1]
        l[i]=float(l[i])
    return l
    fbe.close()    


lbase=beolvaso("output_base_x.csv")
        
def kivono_oszto(xi):
    l=[]
    for i in range(len(xi)):
        l.append((xi[i]-lbase[i])/delta)
    return l    
    
matrix=[]   

for i in range(31):
    fileaktnev="output_x"+str(i)+".csv"
    aktoszlop=beolvaso(fileaktnev)
    modositott_oszlop=kivono_oszto(aktoszlop)
    matrix.append(modositott_oszlop)
    
print(matrix)


fki=open("eredmenymatrix.csv","w")

for i in range(16):
    for j in range(31):
        fki.write(str(matrix[j][i]))
        fki.write(", ")
    fki.write("\n")


fki.close()






