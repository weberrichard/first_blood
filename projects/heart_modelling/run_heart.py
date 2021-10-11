import os

# modell neve
case = "Halasz_P045_heart"

# segéd string a futtatáshoz
s = './heart.out ' + case

# .out kód meghívása
os.system(s)