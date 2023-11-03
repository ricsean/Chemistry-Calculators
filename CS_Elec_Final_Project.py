import FunctionsMod as fn

fn.readFile("Menus.txt", "r", 0, 6) #Display Program Title
while 1:
    try:
        fn.readFile("Menus.txt", "r", 5, 17) #Display Mainmenu
        choice = int(input("\nEnter a Number From the Selection: "))
        if choice == 1: #pH Level
            fn.readFile("Menus.txt", "r", 16, 21) #Display Topic pH Level
            while 1:
                try:
                    fn.readFile("Menus.txt", "r", 20, 33) #Display Menu 
                    ch = int(input("\nEnter a Number From the Selection: "))
                    if ch == 1: #Moles of a Compound
                        fn.readFile("Menus.txt", "r", 32, 43) #Display Title Moles of a Compound
                        ac = input("\nEnter the Acid/Base Compound(Ex. H2O): ")
                        fn.moles(ac)
                    elif ch == 2: #Molarity of an Acid/Base
                        fn.readFile("Menus.txt", "r", 42, 51) #Display Title Molarity of an Acid/Base
                        mol = float(input("\nEnter the Moles of the Compound(mol): "))
                        vol = float(input("Enter the Volume of the Solution(g): "))
                        fn.molarity(mol, vol)
                    elif ch == 3: #Concentration of Hydrogen Ions of an Acid
                        fn.readFile("Menus.txt", "r", 50, 59) #Display Title Concentration of Hydrogen Ions of an Acid
                        M = float(input("\nEnter the Molarity of the Acid(M): "))
                        vol = float(input("Enter the Volume of the Acid(L): "))
                        fn.concentration(M, vol)
                    elif ch == 4: #pH of an Acid Compound
                        fn.readFile("Menus.txt", "r", 58, 73) #Display Title pH of an Acid Compound
                        H = float(input("\nEnter the Concentration of Hydrogen Ions of the Solution(H+): "))
                        fn.pH(H)
                    elif ch == 5: #View History
                        fn.readFileOnly("pHLevelHistory.txt")
                    elif ch == 6: #Delete History
                        fn.deleteFile("pHLevelHistory.txt")
                    elif ch == 0: #Back to Main Menu
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nBack to Main Menu")
                        break
                    else:
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nPlease Enter a Number from the Selection Only!")
                        continue
                except ValueError:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nInvalid Input!")
                    continue
                except:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nSomething Went Wrong")
                    continue
        elif choice == 2: #Strong Acid Titration
            fn.readFile("Menus.txt", "r", 86, 97) #Display Topic Strong Acid Titration
            while 1:
                try:
                    fn.readFile("Menus.txt", "r", 96, 109) #Display Menu 
                    ch = int(input("\nEnter a Number From the Selection: "))
                    if ch == 1: #Acid Volume
                        fn.readFile("Menus.txt", "r", 108, 117) #Display Title Acid Volume
                        acidCon = float(input("\nEnter the Acid Concentration(M): "))
                        baseCon = float(input("Enter the Base Concentration(M): "))
                        baseVol = float(input("Enter the Base Volume(mL): "))
                        fn.findVa(baseVol, baseCon, acidCon)
                    elif ch == 2: #Base Volume
                        fn.readFile("Menus.txt", "r", 116, 125) #Display Title Base Volume
                        acidCon = float(input("\nEnter the Acid Concentration(M): "))
                        baseCon = float(input("Enter the Base Concentration(M): "))
                        acidVol = float(input("Enter the Acid Volume(mL): "))
                        fn.findVb(acidVol, acidCon, baseCon)
                    elif ch == 3: #Base Concentration
                        fn.readFile("Menus.txt", "r", 124, 133) #Display Title Base Concentration
                        acidCon = float(input("\nEnter the Acid Concentration(M): "))
                        baseVol = float(input("Enter the Base Volume(mL): "))
                        acidVol = float(input("Enter the Acid Volume(mL): "))
                        fn.findBc(acidVol, acidCon, baseVol)
                    elif ch == 4: #Acid Concentration
                        fn.readFile("Menus.txt", "r", 132, 141) #Display Title Acid Concentration
                        baseCon = float(input("\nEnter the Base Concentration(M): "))
                        baseVol = float(input("Enter the Base Volume(mL): "))
                        acidVol = float(input("Enter the Acid Volume(mL): "))
                        fn.findAc(baseCon, baseVol, acidVol)
                    elif ch == 5: #View History
                        fn.readFileOnly("StrongAcidTitrationHistory.txt")
                    elif ch == 6: #Delete History
                        fn.deleteFile("StrongAcidTitrationHistory.txt")
                    elif ch == 0: #Back to Main Menu
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nBack to Main Menu")
                        break
                    else:
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nPlease Enter a Number from the Selection Only!")
                        continue
                except ValueError:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nInvalid Input!")
                    continue
                except:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nSomething Went Wrong")
                    continue
        elif choice == 3: #Stoichiometry and Yields
            fn.readFile("Menus.txt", "r", 140, 145) #Display Topic Stoichiometry and Yields
            while 1:
                try:
                    fn.readFile("Menus.txt", "r", 144, 157) #Display Menu 
                    ch = int(input("\nEnter a Number From the Selection: "))
                    if ch == 1: #Balancing Chemical Equation
                        fn.readFile("Menus.txt", "r", 156, 161) #Display Title Balancing Chemical Equation
                        ele1 = input("\nEnter the First Reactant(Ex. H2): ")
                        ele2 = input("Enter the Second Reactant(Ex. O2): ")
                        compound = input("Enter the Product Compound(Ex. H2O): ")
                        fn.balanceComp(ele1, ele2, compound)
                    elif ch == 2: #Theoretical Yield
                        fn.readFile("Menus.txt", "r", 160, 169) #Display Title Theoretical Yield
                        actYield = float(input("\nEnter the Actual Yield(grams): "))
                        perYield = float(input("Enter the Percent Yield(%): "))
                        fn.theoYield(actYield, perYield)
                    elif ch == 3: #Actual Yield
                        fn.readFile("Menus.txt", "r", 168, 177) #Display Title Actual Yield
                        perYield = float(input("\nEnter the Percent Yield(%): "))
                        theoYield = float(input("Enter the Theoretical Yield(grams): "))
                        fn.actualYield(perYield, theoYield)
                    elif ch == 4: #Percent Yield
                        fn.readFile("Menus.txt", "r", 176, 185) #Display Title Percent Yield
                        actYield = float(input("\nEnter the Actual Yield(grams): "))
                        theoYield = float(input("Enter the Theoretical Yield(grams): "))
                        fn.percentYield(actYield, theoYield)
                    elif ch == 5: #View History
                        fn.readFileOnly("StoichiometryYieldHistory.txt")
                    elif ch == 6: #Delete History
                        fn.deleteFile("StoichiometryYieldHistory.txt")
                    elif ch == 0: #Back to Main Menu
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nBack to Main Menu")
                        break
                    else:
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nPlease Enter a Number from the Selection Only!")
                        continue
                except ValueError:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nInvalid Input!")
                    continue
                except:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nSomething Went Wrong")
                    continue
        elif choice == 4: #Ideal and Non-Ideal Gas Laws
            fn.readFile("Menus.txt", "r", 184, 189) #Display Topic Ideal and Non-Ideal Gas Laws
            while 1:
                try:
                    fn.readFile("Menus.txt", "r", 188, 202) #Display Menu 
                    ch = int(input("\nEnter a Number From the Selection: "))
                    if ch == 1: #Ideal Temperature
                        fn.readFile("Menus.txt", "r", 201, 210) #Display Title Ideal Temperature
                        pressure = float(input("\nEnter the Pressure(atm): "))
                        volume = float(input("Enter the Volume(L): "))
                        n = float(input("Enter the Number of Moles of Gas(mol): "))
                        fn.idealTemp(pressure, volume, n)
                    elif ch == 2: #Ideal Pressure
                        fn.readFile("Menus.txt", "r", 209, 218) #Display Title Ideal Pressure
                        n = float(input("\nEnter the Number of Moles of Gas(mol): "))
                        temperature = float(input("Enter the Temperature(K): "))
                        volume = float(input("Enter the Volume(L): "))
                        fn.idealPressure(n, temperature, volume)
                    elif ch == 3: #Ideal Volume
                        fn.readFile("Menus.txt", "r", 217, 226) #Display Title Ideal Volume
                        n = float(input("\nEnter the Number of Moles of Gas(mol): "))
                        temperature = float(input("Enter the Temperature(K): "))
                        pressure = float(input("Enter the Pressure(atm): "))
                        fn.idealVolume(n, temperature, pressure)
                    elif ch == 4: #Non-Ideal Pressure
                        fn.readFile("Menus.txt", "r", 225, 238) #Display Title Non-Ideal Pressure
                        volume = float(input("\nEnter the Volume(L): "))
                        n = float(input("Enter the Number of Moles of Gas(mol): "))
                        a = float(input("Enter the Attraction Between Individual Gas Particles(atm*L^2/mol^2): "))
                        b = float(input("Enter the Average Volume of Individual Gas Particles(L/mol): "))
                        temperature = float(input("Enter the Temperature(K): "))
                        fn.nonIdealPressure(n, temperature, volume, b, a)
                    elif ch == 5: #Difference Between Ideal and Non-Ideal Pressures
                        fn.readFile("Menus.txt", "r", 237, 246) #Display Title Difference Between Ideal and Non-Ideal Conditions
                        ip = float(input("\nEnter the Ideal Pressure(atm): "))
                        nip = float(input("Enter the Non-Ideal Pressure(atm): "))
                        fn.diffIpNip(ip, nip)
                    elif ch == 6: #View History
                        fn.readFileOnly("GasLawsHistory.txt")
                    elif ch == 7: #Delete History
                        fn.deleteFile("GasLawsHistory.txt")
                    elif ch == 0: #Back to Main Menu
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nBack to Main Menu")
                        break
                    else:
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nPlease Enter a Number from the Selection Only!")
                        continue
                except ValueError:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nInvalid Input!")
                    continue
                except:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nSomething Went Wrong")
                    continue
        elif choice == 5: #Conversion of SI Units
            fn.readFile("Menus.txt", "r", 245, 250) #Display Topic Conversion of SI Units
            while 1:
                try:
                    fn.readFile("Menus.txt", "r", 249, 262) #Display Menu 
                    ch = int(input("\nEnter a Number From the Selection: "))
                    if ch == 1: #Mass
                        fn.readFile("Menus.txt", "r", 261, 266) #Display Title Mass Selection
                        fn.readFile("Menus.txt", "r", 265, 276) #Display Selection
                        c = int(input("\nEnter a Number From the Selection: "))
                        if c == 1: #From Milligrams
                            mg = float(input("\nEnter the Millirams(mg): "))
                            fn.milligrams(mg)
                        elif c == 2: #From Grams
                            gram = float(input("\nEnter the Grams(g): "))
                            fn.grams(gram)
                        elif c == 3: #From Kilograms
                            kilogram = float(input("\nEnter the Kilograms(kg): "))
                            fn.kilograms(kilogram)
                        elif c == 4: #From Pounds
                            pound = float(input("\nEnter the Pounds(lb): "))
                            fn.pounds(pound)
                        elif c == 5: #From Ounce
                            ounce = float(input("\nEnter the Ounce(oz): "))
                            fn.ounces(ounce)
                    elif ch == 2: #Pressure
                        fn.readFile("Menus.txt", "r", 275, 280) #Display Title Pressure Selection
                        fn.readFile("Menus.txt", "r", 279, 288) #Display Selection
                        c = int(input("\nEnter a Number From the Selection: "))
                        if c == 1: #From Pascals
                            pa = float(input("\nEnter the Pascals(Pa): "))
                            fn.pascal(pa)
                        elif c == 2: #From Atmospheres
                            atm = float(input("\nEnter the Atmospheres(atm): "))
                            fn.atmosphere(atm)
                        elif c == 3: #From Pound Per Square Inch
                            psi = float(input("\nEnter the Pounds Per Square Inch(psi): "))
                            fn.PSI(psi)
                    elif ch == 3: #Volume
                        fn.readFile("Menus.txt", "r", 287, 292) #Display Title Volume Selection
                        fn.readFile("Menus.txt", "r", 291, 301) #Display Selection
                        c = int(input("\nEnter a Number From the Selection: "))
                        if c == 1: #From Milliliter
                            ml = float(input("\nEnter the Milliliter(ml): "))
                            fn.milliliters(ml)
                        elif c == 2: #From Liter
                            liter = float(input("\nEnter the Liters(L): "))
                            fn.liters(liter)
                        elif c == 3: #From US Gallon
                            gallon = float(input("\nEnter the US Gallon(gal): "))
                            fn.gallons(gallon)
                        elif c == 4: #From Cubic Meters
                            cm = float(input("\nEnter the Cubic Meters(m3): "))
                            fn.cubicmeter(cm)
                    elif ch == 4: #From Temperature
                        fn.readFile("Menus.txt", "r", 300, 305) #Display Title Temperature
                        fn.readFile("Menus.txt", "r", 304, 313) #Display Selection
                        c = int(input("\nEnter a Number From the Selection: "))
                        if c == 1: #From Kelvin
                            k = float(input("\nEnter the Kelvin(K): "))
                            fn.kelvin(k)
                        elif c == 2: #From Celsius
                            C = float(input("\nEnter the Celsius(C): "))
                            fn.celsius(C)
                        elif c == 3: #From Fahrenheit
                            f = float(input("\nEnter the Fahrenheit(F): "))
                            fn.fahrenheit(f)
                    elif ch == 5: #View History
                        fn.readFileOnly("ConversionsHistory.txt")
                    elif ch == 6: #Delete History
                        fn.deleteFile("ConversionsHistory.txt")
                    elif ch == 0: #Back to Main Menu
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nBack to Main Menu")
                        break
                    else:
                        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                        print("\nPlease Enter a Number from the Selection Only!")
                        continue
                except ValueError:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nInvalid Input!")
                    continue
                except:
                    fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
                    print("\nSomething Went Wrong")
                    continue
        elif choice == 0: #Exit Program
            fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
            print("\nThank You for Using the Chemistry Calculators!")
            fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
            break
        else:
            fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
            print("\nPlease Enter a Number from the Selection Only!")
            continue
    except ValueError:
        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
        print("\nInvalid Input!")
        continue
    except:
        fn.readFile("Menus.txt", "r", 0, 1) #Display Division Line
        print("\nSomething Went Wrong")
        continue
