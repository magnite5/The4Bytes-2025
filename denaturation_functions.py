chain_length = df["length"].sum()
molecular_weight = df["molecular_weight"].sum()

count_inst = 0
for i in range(df.shape[0]):
  current_inst = df.loc[i]["instability_index"]
  current_length = df.loc[i]["length"]
  count_inst += current_inst * current_length
Weigthed_instability_index = count_inst/chain_length


count_aromaticity = 0
for i in range(df.shape[0]):
  current_aromaticity = df.loc[i]["aromaticity"]
  current_length = df.loc[i]["length"]
  count_aromaticity += current_aromaticity * current_length
Weigthed_aromaticity = count_aromaticity/chain_length

count_pI = 0
for i in range(df.shape[0]):
  current_pI = df.loc[i]["pI"]
  current_length = df.loc[i]["length"]
  count_pI += current_pI * current_length
Weigthed_pI = count_pI/chain_length


def analyse():
 Strong = 0
 Moderate = 0
 Weak = 0

 value = analyse_chain_length()
 if value == "Strong":
   Strong += 1
 elif value == "Moderate":
   Moderate += 1
 else:
   Weak += 1

 value = analyse_molecular_weight()
 if value == "Strong":
   Strong += 1
 elif value == "Moderate":
    Moderate += 1
 else:
   Weak += 1
 
 value = analyse_aromaticity()
 if value == "Strong":
   Strong += 1
 elif value == "Moderate":
    Moderate += 1
 else:
   Weak += 1
 
 value = analyse_instability_index()
 if value == "Strong":
     Strong += 1
 elif value == "Moderate":
     Moderate += 1
 else:
   Weak += 1
 
 value = analyse_pI()
 if value == "Strong":
     Strong += 1
 elif value == "Moderate":
     Moderate += 1
 else:
   Weak += 1


 if Strong > Moderate and Strong > Weak:
   return "Strong"
 elif Moderate > Strong and Moderate > Weak:
   return "Moderate"
 elif Weak > Strong and Weak > Moderate:
   return "Weak"
 elif Strong == Moderate and Strong > Weak:
   return "Strong"
 elif Strong == Weak and Strong > Moderate:
   return "Strong"
 elif Strong == Weak and Strong == Moderate:
  return "Strong"
 elif Moderate == Weak and Moderate > Strong:
   return "Moderate"
 else:
  return "Strong"


def analyse_chain_length():
  if chain_length > 400:
    return "Strong"
  elif chain_length >= 150 and chain_length <= 400:
    return "Moderate"
  else:
    return "Weak"

def analyse_molecular_weight():
  if molecular_weight > 60000:
    return "Strong"
  elif molecular_weight >= 20000 and molecular_weight <= 60000:
    return "Moderate"
  else:
    return "Weak"

def analyse_aromaticity():
  if Weigthed_aromaticity > 0.10:
    return "Strong"
  elif Weigthed_aromaticity >= 0.07 and Weigthed_aromaticity <= 0.10:
    return "Moderate"
  else:
    return "Weak"

def analyse_instability_index():
  if Weigthed_instability_index < 30:
    return "Strong"
  elif Weigthed_instability_index >= 30 and Weigthed_instability_index <= 40:
    return "Moderate"
  else:
    return "Weak"

def analyse_pI():
  if abs((7- Weigthed_pI)) > 2:
    return "Strong"
  elif abs((7-Weigthed_pI)) >= 1 and abs((7-Weigthed_pI)) <= 2:
    return "Moderate"
  else:
    return "Weak"

print(analyse())
print(analyse_chain_length())
print(analyse_molecular_weight())
print(analyse_aromaticity())
print(analyse_instability_index())
print(analyse_pI())
