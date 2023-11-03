import math as m
import chemlib as cl
import os
from pHcalc.pHcalc import Acid, Neutral, System
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

def readFile(filePath, mode, start, end): #For the Displays Only
    lines = []
    f = open(filePath, mode)
    for i, line in enumerate(f):
        if i in range(start, end):
            lines.append(line.strip())
        elif i > end:
            break
    print()
    for x in lines:
        print(x)
    f.close()

def appendFile(fileName, content): #Create History and Add Contents
    fh = open(fileName, "a")
    fh.write(f"\n{content}\n\n*****************************************************************************************")
    fh.close()

def readFileOnly(fileName): #Reading Histories
    if fileExists(fileName) == True:
        readFile("Menus.txt", "r", 73, 75) #Display Start of History
        f = open(fileName, "r")
        print(f.read())
        f.close()
        readFile("Menus.txt", "r", 75, 76) #Display End of History
    else:
        print("\nNo History Yet!")

def fileExists(fileName): #Checks if a File Exists
    try:
        open(fileName, "r")
        return True
    except:
        return False

def deleteFile(fileName): #Deleting Histories
    readFile("Menus.txt", "r", 76, 83) #Display Delete History Banner
    if fileExists(fileName) == True:
        ask = input("\nRequired Confirmation: ")
        if ask == "OK":
            readFile("Menus.txt", "r", 82, 87) #Display History Deleted
            os.remove(fileName)
        else:
            readFile("Menus.txt", "r", 0, 1) #Display Division Line
            print("\nHistory Deletion Cancelled.")
    else:
        print("\nNo History Yet!")

def divide(x1, x2): #Division
    r = x1/x2
    return r

def moles(compound): #Moles of a Compound
    key = []
    val = []
    AM = []
    I = [] # 'i' is the Molar Mass Output
    mol = cl.Compound(compound).occurences #To Separate the Compound to single elements
    i = 0 
    for keys, vals in mol.items():
        key.append(keys)
        val.append(vals)
        amass = cl.Element(keys).AtomicMass #To Get the Atomic Mass (unique for every element)
        AM.append(amass)
        i = i + (amass * vals)
        I.append(i)
    molarMass = cl.Compound(compound).molar_mass() #Library Shortcut for getting Molar Mass
    if len(mol) == 1: #1 Element
        ln1 = f"\nAnswer:\n\n{compound} = {key[0]}{val[0]}"
        ln2 = f"\nMoles of the Acid/Base Compound(mol) = {AM[0]:.3f} x {val[0]}"
        ln3 = f"Moles of the Acid/Base Compound(mol) = {i:.3f}"
        ln4 = f"The Moles of the Acid/Base Compound is: {molarMass:.3f} mol"
        print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
        appendFile("pHLevelHistory.txt", f"\nMoles of a Compound\n\nInput:\nAcid/Base Compound: {compound}\n\nOutput:\n{ln4}")
    elif len(mol) == 2: #2 Element
        ln1 = f"\nAnswer:\n\n{compound} = {key[0]}{val[0]} + {key[1]}{val[1]}"
        ln2 = f"\nMoles of the Acid/Base Compound(mol) = ({AM[0]:.3f} x {val[0]}) + ({AM[1]:.3f} x {val[1]})"
        ln3 = f"Moles of the Acid/Base Compound(mol) = {I[0]:.3f} + {I[1]:.3f}"
        ln4 = f"Moles of the Acid/Base Compound(mol) = {i:.3f}"
        ln5 = f"The Moles of the Acid/Base Compound is: {molarMass:.3f} mol"
        print(f"{ln1}\n{ln2}\n{ln3}\n{ln4}\n\n{ln5}")
        appendFile("pHLevelHistory.txt", f"\nMoles of a Compound\n\nInput:\nAcid/Base Compound: {compound}\n\nOutput:\n{ln5}")
    elif len(mol) == 3: #3 Element
        ln1 = f"\nAnswer:\n\n{compound} = {key[0]}{val[0]} + {key[1]}{val[1]} + {key[2]}{val[2]}"
        ln2 = f"\nMoles of the Acid/Base Compound(mol) = ({AM[0]:.3f} x {val[0]}) + ({AM[1]:.3f} x {val[1]}) + ({AM[2]:.3f} x {val[2]})"
        ln3 = f"Moles of the Acid/Base Compound(mol) = {I[0]:.3f} + {I[1] - I[0]:.3f} + {I[2] - I[1]:.3f}"
        ln4 = f"Moles of the Acid/Base Compound(mol) = {i:.3f}"
        ln5 = f"The Moles of the Acid/Base Compound is: {molarMass:.3f} mol"
        print(f"{ln1}\n{ln2}\n{ln3}\n{ln4}\n\n{ln5}")
        appendFile("pHLevelHistory.txt", f"\nMoles of a Compound\n\nInput:\nAcid/Base Compound: {compound}\n\nOutput:\n{ln5}")
    elif len(mol) == 4: #4 Element
        ln1 = f"\nAnswer:\n\n{compound} = {key[0]}{val[0]} + {key[1]}{val[1]} + {key[2]}{val[2]} + {key[3]}{val[3]}"
        ln2 = f"\nMoles of the Acid/Base Compound(mol) = ({AM[0]:.3f} x {val[0]}) + ({AM[1]:.3f} x {val[1]}) + ({AM[2]:.3f} x {val[2]}) + ({AM[3]:.3f} x {val[3]})"
        ln3 = f"Moles of the Acid/Base Compound(mol) = {I[0]:.3f} + {I[1] - I[0]:.3f} + {I[2] - I[1]:.3f} + {I[3] - I[2]:.3f}"
        ln4 = f"Moles of the Acid/Base Compound(mol) = {i:.3f}"
        ln5 = f"The Moles of the Acid/Base Compound is: {molarMass:.3f} mol"
        print(f"{ln1}\n{ln2}\n{ln3}\n{ln4}\n\n{ln5}")
        appendFile("pHLevelHistory.txt", f"\nMoles of a Compound\n\nInput:\nAcid/Base Compound: {compound}\n\nOutput:\n{ln5}")
    elif len(mol) == 5: #5 Element
        ln1 = f"\nAnswer:\n\n{compound} = {key[0]}{val[0]} + {key[1]}{val[1]} + {key[2]}{val[2]} + {key[3]}{val[3]} + {key[4]}{val[4]}"
        ln2 = f"\nMoles of the Acid/Base Compound(mol) = ({AM[0]:.3f} x {val[0]}) + ({AM[1]:.3f} x {val[1]}) + ({AM[2]:.3f} x {val[2]}) + ({AM[3]:.3f} x {val[3]}) + ({AM[4]:.3f} x {val[4]})"
        ln3 = f"Moles of the Acid/Base Compound(mol) = {I[0]:.3f} + {I[1] - I[0]:.3f} + {I[2] - I[1]:.3f} + {I[3] - I[2]:.3f} + {I[4] - I[3]:.3f}"
        ln4 = f"Moles of the Acid/Base Compound(mol) = {i:.3f}"
        ln5 = f"The Moles of the Acid/Base Compound is: {molarMass:.3f} mol"
        print(f"{ln1}\n{ln2}\n{ln3}\n{ln4}\n\n{ln5}")
        appendFile("pHLevelHistory.txt", f"\nMoles of a Compound\n\nInput:\nAcid/Base Compound: {compound}\n\nOutput:\n{ln5}")

def molarity(mol, vol): #Molarity of an Acid
    M = divide(vol, mol)
    ln1 = f"\nAnswer:\n\nMolarity(M) = {vol} / {mol}"
    ln2 = f"Molarity(M) = {M:.3f}"
    ln3 = f"The Molarity of the Acid is: {M:.3f} M"
    print(f"{ln1}\n{ln2}\n\n{ln3}")
    appendFile("pHLevelHistory.txt", f"\nMolarity of an Acid\n\nInput:\nMoles of the Compound: {mol} mol\nVolume of the Solution: {vol} g\n\nOutput:\n{ln3}")

def concentration(M, volume): #Concentration of Hydrogen Ions
    c = divide(M, volume)
    ln1 = f"\nAnswer:\n\nConcentration of Hydrogen Ions(H+) = {M} / {volume}"
    ln2 = f"Concentration of Hydrogen Ions(H+) = {c:.3f}"
    ln3 = f"The Concentration of Hydrogen Ions of the Acid is: {c:.3f} H+"
    print(f"{ln1}\n{ln2}\n{ln3}")
    appendFile("pHLevelHistory.txt", f"\nConcentration of Hydrogen Ions of an Acid\n\nInput:\nMolarity of the Acid: {M} M\nVolume of the Solution: {volume} L\n\nOutput:\n{ln3}")

def pH(H): #pH Level of an Acid
    p = abs(m.log10(H)) #pH Level Equation
    if p == 7:
        ln1 = f"\nAnswer:\n\npH = -log({H})"
        ln2 = f"pH = {p:.3f}"
        ln3 = f"The pH Level is: {p:.3f} Nuetral"
        print(f"{ln1}\n{ln2}\n\n{ln3}")
        appendFile("pHLevelHistory.txt", f"\npH of an Acid Compound\n\nInput:\nConcentration of Hydrogen Ions of the Solution: {H} H+\n\nOutput:\n{ln3}")
    elif p > 7:
        ln1 = f"\nAnswer:\n\npH = -log({H})"
        ln2 = f"pH = {p:.3f}"
        ln3 = f"The pH Level is: {p:.3f} Alcaline"
        print(f"{ln1}\n{ln2}\n\n{ln3}")
        appendFile("pHLevelHistory.txt", f"\npH of an Acid Compound\n\nInput:\nConcentration of Hydrogen Ions of the Solution: {H} H+\n\nOutput:\n{ln3}")
    elif p < 7:
        ln1 = f"\nAnswer:\n\npH = -log({H})"
        ln2 = f"pH = {p:.3f}"
        ln3 = f"The pH Level is: {p:.3f} Acidic"
        print(f"{ln1}\n{ln2}\n\n{ln3}")
        appendFile("pHLevelHistory.txt", f"\npH of an Acid Compound\n\nInput:\nConcentration of Hydrogen Ions of the Solution: {H} H+\n\nOutput:\n{ln3}")

def titrations(a, b, c): #Titration Equation
    r = (a * b) / c
    return r

def findVa(baseVol, baseCon, acidCon): #Find Acid Volume
    Va = titrations(baseVol, baseCon, acidCon)
    ln1 = f"\nAnswer:\n\nAcid Volume(Va) = ({baseVol} x {baseCon}) / {acidCon}"
    ln2 = f"Acid Volume(Va) = {baseVol * baseCon:.3f} / {acidCon}"
    ln3 = f"Acid Volume(Va) = {Va:.3f}"
    ln4 = f"The Acid Volume is: {Va:.3f} mL"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StrongAcidTitrationHistory.txt", f"\nAcid Volume\n\nInput:\nAcid Concentration[H+]: {acidCon} M\nBase Concentration[OH-]: {baseCon} M\nBase Volume: {baseVol} mL\n\nOutput:\n{ln4}")
    curve(baseCon, baseVol, Va)

def findVb(acidVol, acidCon, baseCon): # Find Base Colume
    Vb = titrations(acidVol, acidCon, baseCon)
    ln1 = f"\nAnswer:\n\nBase Volume(Vb) = ({acidVol} x {acidCon}) / {baseCon}"
    ln2 = f"Base Volume(Vb) = {acidVol * acidCon:.3f} / {baseCon}"
    ln3 = f"Base Volume(Vb) = {Vb:.3f}"
    ln4 = f"The Base Volume is: {Vb:.3f} mL"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StrongAcidTitrationHistory.txt", f"\nBase Volume\n\nInput:\nAcid Concentration[H+]: {acidCon} M\nBase Concentration[OH-]: {baseCon} M\nAcid Volume: {acidVol} mL\n\nOutput:\n{ln4}")
    curve(baseCon, Vb, acidVol)

def findBc(acidVol, acidCon, baseVol): #Find Base Concentration
    Bc = titrations(acidVol, acidCon, baseVol)
    ln1 = f"\nAnswer:\n\nBase Concentration[OH-] = ({acidVol} x {acidCon}) / {baseVol}"
    ln2 = f"Base Concentration[OH-] = {acidVol * acidCon:.3f} / {baseVol}"
    ln3 = f"Base Concentration[OH-] = {Bc:.3f}"
    ln4 = f"The Base Concentration[OH-] is: {Bc:.3f} M"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StrongAcidTitrationHistory.txt", f"\nBase Concentration[OH-]\n\nInput:\nAcid Concentration[H+]: {acidCon} M\nBase Volume: {baseVol} mL\nAcid Volume: {acidVol} mL\n\nOutput:\n{ln4}")
    curve(Bc, baseVol, acidVol)

def findAc(baseVol, baseCon, acidVol): #Find Base Concentration
    Ac = titrations(baseVol, baseCon, acidVol)
    ln1 = f"\nAnswer:\n\nAcid Concentration[H+] = ({baseVol} x {baseCon}) / {acidVol}"
    ln2 = f"Acid Concentration[H+] = {baseVol * baseCon:.3f} / {acidVol}"
    ln3 = f"Acid Concentration[H+] = {Ac:.3f}"
    ln4 = f"The Acid Concentration[H+] is: {Ac:.3f} M"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StrongAcidTitrationHistory.txt", f"\nAcid Concentration[H+]\n\nInput:\nBase Concentration[OH-]: {baseCon} M\nBase Volume: {baseVol} mL\nAcid Volume: {acidVol} mL\n\nOutput:\n{ln4}")
    curve(baseCon, baseVol, acidVol)

def curve(H, bVol, aVol):
    na_moles = np.linspace(1e-8, H, 50)
    sol_volume = bVol + aVol / 1000 # Liter
    phos = Acid(pKa=[H], charge=1, conc=1.e-3)
    phs = []
    for mol in na_moles:
        na = Neutral(charge=1, conc=mol/sol_volume)
        system = System(phos, na)
        system.pHsolve(guess_est=True)
        phs.append(system.pH)
    plt.plot(na_moles, phs, marker = "*")
    plt.xlabel("Moles")
    plt.ylabel("pH Level")
    plt.show()

def balanceComp(ele1, ele2, compound): #Balancing Chemical Equation
    r = cl.Reaction(reactants = [cl.Compound(ele1), cl.Compound(ele2)], products = [cl.Compound(compound)]) #to Compile the Reactants and Product
    if r.is_balanced == True: #To Check if Balanced
        ln1 = "\nAnswer:\n\nThe Entered Chemical Stoichiometry Equation:"
        ln2 = r.formula
        ln3 = "is Already Balanced."
        ln4 = "Note: The Stoichiometric Coefficient is the Number Before the \nChemical Formula in a Balanced Equation."
        lines = f"\n{ln1}\n{ln2}\n{ln3}\n\n{ln4}"
        print(lines)
        appendFile("StoichiometryYieldHistory.txt", f"\nBalancing Chemical Equation\n{lines}")
    else:
        ln1 = "\nAnswer:\n\nThe Entered Chemical Stoichiometry Equation:"
        ln2 = r.formula
        ln3 = "is Not Balanced."
        ln4 = "The Balanced Chemical Stoichiometry Equation is:"
        r.balance() #To Balance the Chemical Equation
        ln5 = r.formula
        ln6 = "Note: The Stoichiometric Coefficient is the Number Before the \nChemical Formula in a Balanced Equation."
        lines = f"\n{ln1}\n{ln2}\n{ln3}\n\n{ln4}\n{ln5}\n\n{ln6}"
        print(lines)
        appendFile("StoichiometryYieldHistory.txt", f"\nBalancing Chemical Equation\n{lines}")

def yields(x, y, z): #Yields Equation
    y = (x / y) * z
    return y

def theoYield(actual, percent): #Find Theoretical Yield
    theo = yields(actual, percent, 100)
    ln1 = f"\nAnswer:\n\nTheoretical Yield: ({actual} / {percent}) x 100"
    ln2 = f"Theoretical Yield: {actual / percent:.3f} x 100"
    ln3 = f"Theoretical Yield: {theo:.3f} g"
    ln4 = f"The Theoretical Yield is: {theo:.3f} g"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StoichiometryYieldHistory.txt", f"Theoretical Yield\n\nInput:\nActual Yield: {actual} g\nPercent Yield: {percent}%\n\nOutput:\n{ln4}")

def actualYield(percent, theoretical): #Find Actual Yield
    act = yields(percent, 100, theoretical)
    ln1 = f"\nAnswer:\n\nActual Yield: ({percent} / 100) x {theoretical}"
    ln2 = f"Actual Yield: {percent / 100:.3f} x {theoretical}"
    ln3 = f"Actual Yield: {act:.3f} g"
    ln4 = f"The Actual Yield is: {act:.3f} g"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StoichiometryYieldHistory.txt", f"Actual Yield\n\nInput:\nPercent Yield: {percent}%\nTheoretical Yield: {theoretical} g\n\nOutput:\n{ln4}")

def percentYield(actual, theoretical): #Find Percent Yield
    per = yields(actual, theoretical, 100)
    ln1 = f"\nAnswer:\n\nPercent Yield: ({actual} / {theoretical}) x 100"
    ln2 = f"Percent Yield: {actual / theoretical:.3f} x 100"
    ln3 = f"Percent Yield: {per:.3f} g"
    ln4 = f"The Percent Yield is: {per:.3f}%"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("StoichiometryYieldHistory.txt", f"Percent Yield\n\nInput:\nActual Yield: {actual} g\nTheoretical Yield: {theoretical} g\n\nOutput:\n{ln4}")

def idealTemp(pressure, volume, n): #Find Ideal Temperature
    R = 8.20578 * m.pow(10, -2) #Ideal Gas Law Constant
    it = (pressure * volume) / (n * R)
    ln1 = f"\nAnswer:\n\nIdeal Temperature: ({pressure} atm x {volume} L) / ({n} mol x {R:.3f} L atm K-1 mol-1)"
    ln2 = f"Ideal Temperature: {pressure * volume:.3f} / ({n * R:.3f})"
    ln3 = f"Ideal Temperature: {it:.3f} K"
    ln4 = f"The Ideal Temperature is: {it:.3f} K"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("GasLawsHistory.txt", f"Ideal Temperature\n\nInput:\nPressure: {pressure} atm\nVolume: {volume} L\nNumber of Moles of Gas: {n} mol\n\nOutput:\n{ln4}")

def idealPressure(n, temp, vol): #Find Ideal Pressure
    R = 8.20578 * m.pow(10, -2) #Ideal Gas Law Constant
    ip = (n * R * temp) / vol #Equation
    ln1 = f"\nAnswer:\n\nIdeal Pressure: ({n} mol x {R:.3f} L atm K-1 mol-1 x {temp}) / {vol} L"
    ln2 = f"Ideal Pressure: {n * R * temp:.3f}) / {vol} L"
    ln3 = f"Ideal Pressure: {ip:.3f} atm"
    ln4 = f"The Ideal Pressure is: {ip:.3f} atm"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("GasLawsHistory.txt", f"Ideal Pressure\n\nInput:\nNumber of Moles of Gas: {n} mol\nTemperature: {temp} K\nVolume: {vol} L\n\nOutput:\n{ln4}")

def idealVolume(n, temp, pres): #Find Ideal Volume
    R = 0.08206 #Ideal Gas Law Constant
    ip = (n * R * temp) / pres #Equation
    ln1 = f"\nAnswer:\n\nIdeal Volume: ({n} mol x {R:.3f} L atm K-1 mol-1 x {temp}) / {pres} atm"
    ln2 = f"Ideal Volume: {n * R * temp:.3f}) / {pres} atm"
    ln3 = f"Ideal Volume: {ip:.3f} L"
    ln4 = f"The Ideal Volume is: {ip:.3f}  L"
    print(f"{ln1}\n{ln2}\n{ln3}\n\n{ln4}")
    appendFile("GasLawsHistory.txt", f"Ideal Volume\n\nInput:\nNumber of Moles of Gas: {n} mol\nTemperature: {temp} K\nPressure: {pres} L\n\nOutput:\n{ln4}")

def nonIdealPressure(n, temp, vol, b, a): #Find Non-Ideal Pressure
    R = 8.20578 * m.pow(10, -2) #Ideal Gas Law Constant
    x = (n * R * temp) / (vol - (n * b)) #Non Ideal Pressure Equation Divided into 2 parts x and y
    y = a * m.pow((n/vol), 2)
    nip = x - y
    ln1 = f"\nAnswer:\n\nNon-Ideal Pressure: X = \n(({n} mol x {R:.3f} L*atm/mol*K x {temp} K) / ({vol} L - ({n} mol x {b} L/mol)))"
    ln2 = f"Non-Ideal Pressure: X = ({n * R * temp:.3f} / ({vol} L - {n * b:.3f}))"
    ln3 = f"Non-Ideal Pressure: X = {n * R * temp:.3f} / {vol - (n * b):.3f}"
    ln4 = f"Non-Ideal Pressure: X = {x:.3f}"
    ln5 = f"Non-Ideal Pressure: Y = {a} atm*L^2/mol^2 x ({n} mol / {vol} L)^2"
    ln6 = f"Non-Ideal Pressure: Y = {a} atm*L^2/mol^2 x ({n / vol:.3f})^2)"
    ln7 = f"Non-Ideal Pressure: Y = {a} atm*L^2/mol^2 x {m.pow((n * vol), 2):.3f}"
    ln8 = f"Non-Ideal Pressure: Y = {y:.3f})"
    ln9 = f"Non-Ideal Pressure: X - Y = {x:.3f} - {y:.3f}"
    ln10 = f"Non-Ideal Pressure: {nip:.3f} atm"
    ln11 = f"The Non-Ideal Pressure is: {nip:.3f} atm"
    print(f"{ln1}\n{ln2}\n{ln3}\n{ln4}\n\n{ln5}\n{ln6}\n{ln7}\n{ln8}\n{ln9}\n{ln10}\n\n{ln11}")
    appendFile("GasLawsHistory.txt", f"Non-Ideal Pressure\n\nInput:\nVolume: {vol} L\nNumber of Moles of Gas: {n} mol\nAttraction Between Individual Gas Particles: {a} atm*L^2/mol^2\nAverage Volume of Individual Gas Particles: {b} L/mol\nTemperature: {temp} K\n\nOutput:\n{ln11}")

def diffIpNip(ip, nip): #Find Difference Between Ideal and Non-Ideal Pressures
    diff = abs(ip - nip)
    if ip < nip:
        line = f"\nAnswer:\n\nThe pressure for the ideal gas is {ip} atm and the\npressure of the non-ideal gas is {nip} atm.\nThe non-ideal gas have a greater pressure by {diff:.3f} atm."
        print(line)
        appendFile("GasLawsHistory.txt", f"Difference Between Ideal and Non-Ideal Pressures\n{line}")
    elif ip > nip:
        line = f"\nAnswer:\n\nThe pressure for the ideal gas is {ip} atm and the\npressure of the non-ideal gas is {nip} atm.\nThe ideal gas have a greater pressure by {diff:.3f} atm."
        print(line)
        appendFile("GasLawsHistory.txt", f"Difference Between Ideal and Non-Ideal Pressures\n{line}")
    elif ip == nip:
        line = f"\nAnswer:\n\nThe pressure for the ideal gas is {ip} atm and the\npressure of the non-ideal gas is {nip} atm.\nBoth have the same pressure."
        print(line)
        appendFile("GasLawsHistory.txt", f"Difference Between Ideal and Non-Ideal Pressures\n{line}")

def conv(num, cons1, cons2): #Most Conversions Formula
    output = (num * cons1) / cons2
    return output

def milligrams(mg): #Convert Milligrams to other Metrics of Mass
    mgTog = mg / 1000
    mgTokg = mg / 1000000
    mgTolb = mg / 453592.37
    mgTooz = mg / 28349.5231
    ln = f"\nMilligrams to Grams: {mgTog}g\nMilligrams To Kilograms: {mgTokg}kg\nMilligrams to Pounds: {mgTolb}lb\nMilligrams to Ounce: {mgTooz}oz"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Mass Milligram Conversions\n{ln}")

def grams(g): #Convert Grams to other Metrics of Mass
    gTomg = g * 1000
    gTokg = g * constants.gram
    gTolb = conv(g, constants.gram, constants.lb)
    gTooz = conv(g, constants.gram, constants.oz)
    ln = f"\nGrams to Milligrams: {gTomg}mg\nGrams To Kilograms: {gTokg}kg\nGrams to Pounds: {gTolb}lb\nGrams to Ounce: {gTooz}oz"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Mass Gram Conversions\n{ln}")

def kilograms(kg): #Convert Kilograms to other Metrics of Mass
    kgTomg = kg * 1000000
    kgTog = kg / constants.gram
    kgTolb = kg / constants.lb
    kgTooz = kg / constants.oz
    ln = f"\nKilograms to Milligrams: {kgTomg}mg\nKilograms To Grams: {kgTog}g\nKilograms to Pounds: {kgTolb}lb\nKilograms to Ounce: {kgTooz}oz"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Mass Kilogram Conversions\n{ln}")

def pounds(lb): #Convert Pounds to other Metrics of Mass
    lbTomg = (lb * constants.lb) * 1000000
    lbTog = conv(lb, constants.lb, constants.gram)
    lbTokg = lb * constants.lb
    lbTooz = conv(lb, constants.lb, constants.oz)
    ln = f"\nPounds to Milligrams: {lbTomg}mg\nPounds To Grams: {lbTog}g\nPounds to Kilograms: {lbTokg}kg\nPounds to Ounce: {lbTooz}oz"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Mass Pound Conversions\n{ln}")

def ounces(oz): #Convert Ounces to other Metrics of Mass
    ozTomg = (oz * constants.oz) * 1000000
    ozTog = conv(oz, constants.oz, constants.gram)
    ozTokg = oz * constants.oz
    ozTolb = conv(oz, constants.oz, constants.lb)
    ln = f"\nOunce to Milligrams: {ozTomg}mg\nOunce To Grams: {ozTog}g\nOunce to Kilograms: {ozTokg}lb\nOunce to Pounds: {ozTolb}lb"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Mass Ounce Conversions\n{ln}")

def pascal(pa): #Convert Pascals to other Metrics of Mass
    paToatm = pa / constants.atm
    paTopsi = pa / constants.psi
    ln = f"\nPascals to Atmosphere: {paToatm}atm\nPascals to Pounds Per Square Inch: {paTopsi}psi"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Pressure Pascals Conversions\n{ln}")

def atmosphere(atm): #Convert Atmospheres to other Metrics of Mass
    atmTopa = atm * constants.atm
    atmTopsi = conv(atm, constants.atm, constants.psi)
    ln = f"\nAtmosphere to Pascals: {atmTopa}pa\nAtmosphere to Pounds Per Square Inch: {atmTopsi}psi"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Pressure Atmosphere Conversions\n{ln}")

def PSI(psi): #Convert Pounds Per Square Inch to other Metrics of Mass
    psiTopa = psi * constants.psi
    psiToatm = conv(psi, constants.psi, constants.atm)
    ln = f"\nPounds Per Square Inch to Pascals: {psiTopa}pa\nPounds Per Square Inch to Atmosphere: {psiToatm}atm"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Pressure Pounds Per Square Inch Conversions\n{ln}")

def milliliters(ml): #Convert Milliliters to other Metrics of Mass
    mlTol = ml / 1000
    mlTogal = ml / 3785.41178
    mlTocm = ml / 1000000
    ln = f"\nMilliliters to Liters: {mlTol}L\nMilliliters To US Gallon: {mlTogal}gal\nMilliliters to Cubic Meters: {mlTocm}m3"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Volume Milliliter Conversions\n{ln}")

def liters(l): #Convert Liter to other Metrics of Mass
    lToml = l * 1000
    lTogal = conv(l, constants.liter, constants.gallon)
    lTocm = l * constants.liter
    ln = f"\nLiters to Milliliters: {lToml}ml\nLiters To US Gallon: {lTogal}gal\nLiters to Cubic Meters: {lTocm}m3"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Volume Liter Conversions\n{ln}")

def gallons(gal): #Convert US Gallon to other Metrics of Mass
    galToml = gal * 3785.41178
    galTol = conv(gal, constants.gallon, constants.liter)
    galTocm = gal * constants.gallon
    ln = f"\nUS Gallon to Milliliters: {galToml}ml\nUS Gallon To Liters: {galTol}L\nUS Gallons to Cubic Meters: {galTocm}m3"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Volume US Gallon Conversions\n{ln}")

def cubicmeter(cm): #Convert Cubic Meter to other Metrics of Mass
    cmToml = cm * 1000000
    cmTol = cm / constants.liter
    cmTogal = cm / constants.gallon
    ln = f"\nCubic Meters to Milliliters: {cmToml}mg\nCubic Meters To Liters: {cmTol}g\nCubic Meters to US Gallon: {cmTogal}gal"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Volume Cubic Meter Conversions\n{ln}")

def kelvin(K): #Convert Kelvin to other Metrics of Mass
    ktoc = K - constants.zero_Celsius
    ktof = (K - 273.15) * 1.8000 + 32.00
    ln = f"\nKelvin to Celsius: {ktoc}C\nKelvin to Fahrenheit: {ktof}F"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Temperature Kelvin Conversions\n{ln}")

def celsius(C): #Convert Celsius to other Metrics of Mass
    cTok = C + constants.zero_Celsius
    cTof = C * 1.8000 + 32.00
    ln = f"\nCelsius to Kelvin: {cTok}K\nCelsius to Fahrenheit: {cTof}F"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Temperature Celsius Conversions\n{ln}")

def fahrenheit(F): #Convert Fahrenheit to other Metrics of Mass
    fTok = (F - 32) / 1.8 + constants.zero_Celsius
    fToc = (F - 32) * constants.degree_Fahrenheit
    ln = f"\nFahrenheit to Kelvin: {fTok}K\nFahrenheit to Celsius: {fToc}C"
    print(f"\nConversions\n{ln}")
    appendFile("ConversionsHistory.txt", f"Temperature Fahrenheit Conversions\n{ln}")