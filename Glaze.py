# Glaze calculator
import numpy as np
import re

class Element:
    def __init__(self, symbol, name, atWt):
        self.symbol = symbol
        self.name = name
        self.atWt = atWt

element_data = [ \
    ('Aluminum', 'Al', 26.97), \
    ('Barium', 'Ba', 137.36), \
    ('Bismuth', 'Bi', 209.00), \
    ('Boron', 'B', 10.82), \
    ('Cadmium', 'Cd', 112.42), \
    ('Calcium', 'Ca', 40.08),\
    ('Carbon', 'C', 12.01),\
    ('Chlorine', 'Cl', 35.457),\
    ('Chromium', 'Cr', 52.01),\
    ('Cobalt', 'Co', 58.94),\
    ('Copper', 'Cu', 63.54),\
    ('Flourine', 'F', 19.00),\
    ('Gold', 'Au', 197.20),\
    ('Hydrogen', 'H', 1.008),\
    ('Iron', 'Fe', 55.84),\
    ('Lead', 'Pb', 207.21),\
    ('Lithium', 'Li', 6.94),\
    ('Magnesium', 'Mg', 24.32),\
    ('Manganese', 'Mn', 54.93),\
    ('Molybdenum', 'Mo', 95.98),\
    ('Nickel', 'Ni', 58.69),\
    ('Nitrogen', 'N', 14.008),\
    ('Oxygen', 'O', 16.00),\
    ('Phosphorus', 'P', 30.98),\
    ('Potassium', 'K', 39.096),\
    ('Silicon', 'Si', 28.06),\
    ('Silver', 'Ag', 107.88),\
    ('Sodium', 'Na', 22.997),\
    ('Sulphur', 'S', 36.006),\
    ('Tin', 'Sn', 118.70),\
    ('Titanium', 'Ti', 47.90),\
    ('Vanadium', 'V', 50.95),\
    ('Zinc', 'Zn', 65.38),\
    ('Zirconium', 'Zr', 91.22)
    ]

mpSymElem = {}

for name, symbol, atwt in element_data:
    elem = Element(symbol, name, atwt)
    mpSymElem[symbol] = elem

# print mpSymElem['Al'].name

class StringTokens:
    "Tokenizer"
    
    def __init__(self, s):
        self.s = s
        self.iCh = 0

    def is_empty(self):
        return self.iCh >= len(self.s)

    def peek(self):
        ch = '\a' # something safe to dereference but not valid
        if not(self.is_empty()):
            ch = self.s[self.iCh]
        #print 'peek =', ch
        return ch

    def pop(self):
        ch = self.peek()
        if not(self.is_empty()):
            self.iCh += 1
        #print 'pop =', ch
        return ch


class Molecule:
    "A single oxide, Eg. 'Fe2O3'"
    
    mpNameMol = {}
    nSort = 0
    
    def __init__(self, name, molWt, nSort):

        self.name = name
        self.molWt = molWt
        self.nSort = nSort

    @classmethod
    def ensure(cls, name):
        if name in Molecule.mpNameMol:
            return Molecule.mpNameMol[name]

        # mpElemC = mols of elem / mol of material
        # mpElemC / molWt = mols of elem / gram material
        mpElemC = {}
        global mpSymElem

        toks = StringTokens(name)

        while not(toks.is_empty()):
            ch = toks.pop()
            if ch.isupper():
                sym = ch
                if toks.peek().islower():
                    sym += toks.pop()
                elem = mpSymElem[sym]
                c = 1
                if toks.peek().isdigit():
                    c = int(toks.pop())
                # print mul * c, " ", elem.name
                mpElemC[elem] = mpElemC.get(elem, 0) + c
            else:
                raise Exception('unrecognized character in formula')
            
        # Calculate molecular weight as sum of atomic weights        
        molWt = 0
        for elem, c in mpElemC.items():
            molWt += c * elem.atWt

        mol = cls(name, molWt, Molecule.nSort)
        Molecule.nSort += 1

        Molecule.mpNameMol[name] = mol
        
        return mol

##    Python 2?
##    def __cmp__(self, other):
##        return self.nSort - other.nSort
    
    def __lt__(self, other):
        return self.nSort < other.nSort

    def __hash__(self):
        return hash(self.nSort)


for name in ['SiO2', 'Al2O3', 'Na2O', 'K2O', 'MgO', 'CaO', 'TiO2', 'Fe2O3']:
    Molecule.ensure(name)


        
class Formula:
    "A mixture of oxides, Eg. 'Al2O3 2SiO2'"
    
    def __init__(self, stForm):
        self.name = stForm
        
        # mpMolC = mols of oxide / mol of material
        # mpMolC / molWt = mols of oxide / gram material
        self.mpMolC = {}
        global mpSymElem

        toks = StringTokens(stForm)

        mul = 1
        while not(toks.is_empty()):
            ch = toks.pop()
            if ch.isdigit():
                mul = float(ch)
                while toks.peek().isdigit():
                    mul = 10 * mul + float(toks.pop())
                if toks.peek() == '.':
                    toks.pop()
                    digit_value = 0.1
                    while toks.peek().isdigit():
                        mul += float(toks.pop()) * digit_value
                        digit_value /= 10
            elif ch == ' ' or ch == '.':
                mul = 1
            elif ch.isupper():
                sym = ch
                while toks.peek().isalnum():
                    sym += toks.pop()
                mol = Molecule.ensure(sym)
                # print mul " ", mol.name
                self.mpMolC[mol] = self.mpMolC.get(mol, 0) + mul
            else:
                raise Exception('unrecognized character in formula')
            
        # Calculate molecular weight as sum of atomic weights        
        molWt = 0
        for mol, c in self.mpMolC.items():
            molWt += c * mol.molWt
        self.molWt = molWt
        

class Material:
    "A raw material for a glaze"
    
    def __init__(self, names, stRaw, stFired):

        if type(names) == str:
            self.names = [names]
        else:
            self.names = names

        self.fAdd = (self.names[0][0] == '+')
        
        # strip +'s
        for i,name in enumerate(self.names):
            if name[0] == '+':
                self.names[i] = name[1:]

        self.name = self.names[0]

        self.mpMolFrac = {}

        if not stFired:
            return
        
        formRaw = Formula(stRaw)

        formFired = Formula(stFired)
        
        for mol, c in formFired.mpMolC.items():
            # How many mols of this oxide per gram of material
            self.mpMolFrac[mol] = float(c) * mol.molWt / formRaw.molWt
        
        self.calcLoi()

    @classmethod
    def fromPercents(cls, names, percents):
        "Initialize Material from an array of formula/percentage pairs"

        # BB can I do this constructor via optional param and avoid calling regular init?
        
        material = cls(names, '', '')

        for sForm, percent in percents:
            if sForm is "LOI":
                continue
            form = Formula(sForm)
            frac = float(percent) / 100.0

            for mol, c in form.mpMolC.items():
                assert(c == 1)
                # How many mols of this oxide per gram of material
                material.mpMolFrac[mol] = frac # handle dups?

        material.calcLoi()
            
        return material

    def calcLoi(self):
        fracFired = 0.0
        for frac in self.mpMolFrac.values():
            fracFired += frac
        
        self.loi = 1.0 - fracFired

        self.mpMolFracFired = {}
        for mol,frac in self.mpMolFrac.items():
            self.mpMolFracFired[mol] = frac / fracFired
        
        #print name, self.mpMolFrac, "loi =", self.loi


material_data = [
    ('Whiting', "CaCO3", "CaO"),
#    ('Kaolin', 'Al2O3 2SiO2 2H2O', 'Al2O3 2SiO2'),
    ('Talc', '3MgO 4SiO2 H2O', '3MgO 4SiO2'),
    (['Silica', "Flint"], 'SiO2', 'SiO2'),
    ('Zinc Oxide', 'ZnO', 'ZnO'),
    ('+RIO', 'Fe2O3', 'Fe2O3'),
    ('+Copper Carbonate', 'CuCO3', 'CuO'),
    ('+Titanium Dioxide', 'TiO2', 'TiO2'),
    (['CuO', 'Black Copper Oxide'], 'CuO', 'CuO'),
    ]

Materials = {}
def addMat(mat):
    for name in mat.names:
        Materials[name] = mat
        Materials["+" + name]= mat
    
for names, stRaw, stFired in material_data:
    addMat(Material(names, stRaw, stFired))
    #print "mol wt of", material.name, "is", material.formRaw.molWt

addMat(Material.fromPercents('Barnard Clay',\
            [('SiO2', 46.99),\
             ('Al2O3', 7.60),\
             ('MgO', 0.70), ('CaO', 0.60), ('K2O', 1.10), ('Na2O', 0.6), ('Fe2O3', 33.9), ('TiO2', 0.2)]))
addMat(Material.fromPercents(['Custer Feldspar', "Custer", "Potash Feldspar"],\
            [('SiO2', 68.5),\
             ('Al2O3', 17.0),\
             ('CaO', 0.30), ('K2O', 10.0), ('Na2O', 3.00), ('Fe2O3', 0.1)]))
addMat(Material.fromPercents(['EPK', "EP Kaolin", "Kaolin"],\
            [('SiO2', 45.73),\
             ('Al2O3', 37.36),\
             ('MgO', 0.10), ('CaO', 0.18), ('K2O', 0.35), ('Na2O', 0.07), ('P2O5', 0.26), ('Fe2O3', 0.76), ('TiO2', 0.38)]))
addMat(Material.fromPercents(['Cornwall Stone', 'Cornish Stone'],\
            [('SiO2', 73.78),
             ('Al2O3', 16.33),
             ('MgO', 0.14), ('CaO', 1.81), ('K2O', 4.30), ('Na2O', 3.30), ('Fe2O3', 0.19), ('TiO2', 0.15)]))
addMat(Material.fromPercents('Dolomite',\
            [('CaO', 30.48), ('MgO', 21.91)]))
addMat(Material.fromPercents(['+Rutile'],\
            [('Fe2O3', 9.90), ('TiO2', 90.0)]))
addMat(Material.fromPercents('Nepheline Syenite',\
            [('SiO2', 60.70),\
             ('Al2O3', 23.30),\
             ('MgO', 0.10), ('CaO', 0.70), ('K2O', 4.60), ('Na2O', 9.80), ('Fe2O3', 0.10)]))
addMat(Material.fromPercents('+Bentonite',\
            [('SiO2', 59.33),\
             ('Al2O3', 20.12),\
             ('MgO', 2.01), ('CaO', 1.01), ('K2O', 1.01), ('Na2O', 3.02), ('Fe2O3', 3.51)]))

addMat(Material.fromPercents("Gerstley Borate",\
            [('SiO2', 14.80),\
             ('Al2O3', 1.00),\
             ('B2O3', 26.80), ('MgO', 3.50), ('CaO', 19.40), ('K2O', 0.40), ('Na2O', 4.00), ('P2O5', 0.10), ('Fe2O3', 0.40), ('TiO2', 0.10)]))

addMat(Material.fromPercents(["Kona F-4", "Kona", "Soda Spar"],\
            [('SiO2', 66.74),\
             ('Al2O3', 19.58),\
             ('MgO', 0.05), ('CaO', 1.70), ('K2O', 4.80), ('Na2O', 6.89), ('Fe2O3', 0.04)]))

addMat(Material.fromPercents(["OM-4 Ball Clay", "OM4", "Ball Clay", "OM-4"],\
            [('CaO', 0.20),\
             ('K2O', 1.10),\
             ('MgO', 0.50),\
             ('Na2O', 0.20),\
             ('TiO2', 2.00),\
             ('Al2O3', 26.30),\
             ('SiO2', 59.70),\
             ('Fe2O3', 1.30)]))

addMat(Material.fromPercents("Petalite",\
            [('Li2O', 4.87),\
             ('Al2O3', 16.65),\
             ('SiO2', 78.49)]))

addMat(Material.fromPercents(["Zircopax", "Zirconium Silicate"],
            [('ZrO2', 67.21),
             ('SiO2', 32.79)]))

addMat(Material.fromPercents(["Wood Ash", "Ash"],
            [('CaO', 37),
             ('MgO', 11),
             ('K2O', 16),
             ('Na2O', 1),
             ('P2O5', 9),
             ('SiO2', 7),
             ('Fe2O3', 1),
             ('MnO', 2),
             ('LOI', 16)]))

addMat(Material.fromPercents(["Yellow Ochre", "Ochre"],
             [("CaO", 0.23),
              ("MgO", 0.06),
              ("Al2O3", 19.92),
              ("SiO2", 57.78),
              ("Fe2O3", 22.01)]))


#print Molecule.mpNameMolWt

##for mat in Materials:
##    print mat, Materials[mat].mpMolFrac

def mulOxide(mpMolWt, sMol, gMul):
    mol = Molecule.ensure(sMol)
    if mol in mpMolWt:
        mpMolWt[mol] *= gMul
    else:
        raise("mulOxide on missing oxide")

def replaceOxide(mpMolWt, sMolOld, frac, sMolNew):
    molOld = Molecule.ensure(sMolOld)
    molNew = Molecule.ensure(sMolNew)

    if molOld in mpMolWt:
        wtOld = mpMolWt[molOld]
        if frac < 1.0:
            mpMolWt[molOld] = wtOld * (1 - frac)
        else:
            del mpMolWt[molOld]
        nOld = wtOld / molOld.molWt
        nNew = nOld * frac
        wtNew = nNew * molNew.molWt
        mpMolWt[molNew] = mpMolWt.get(molNew, 0) + wtNew

class Glaze:
    def __init__(self, ingreds, name = None):
        self.name = name
        self.mats = []      # ingredients in recipe order
        self.mpMatWt = {}   # how much of each raw ingredient in 100g raw glaze
        self.mpMatFAdd = {}

        self.wtTotal = 0.0  # total weight of 100g batch, including additives
        
        # BB should this just be a Formula? Rationalize percent vs. frac
        self.mpMolWt = {} # grams fired oxide/100g raw glaze
        self.mpMolWtFired = {} # grams fired oxide/100g fired glaze

        if not ingreds:
            return
        
        sum = 0.0
        sumAdditives = 0.0
        for sMat, wt in ingreds:
            fAdd = False
            if sMat[0] == "+":
                fAdd = True
                sMat = sMat[1:]
            else:
                sum += wt
            sumAdditives += wt
            self.wtTotal += wt
            mat = Materials[sMat]
            self.mats.append(mat)
            self.mpMatWt[mat] = wt
            self.mpMatFAdd[mat] = fAdd

        # normalize to 100

        if sum == 0.0:
            sum = sumAdditives
            
        for mat in self.mpMatWt:
            self.mpMatWt[mat] *= 100.0 / sum
        self.wtTotal *= 100.0 / sum
            
        for mat,wt in self.mpMatWt.items():
            for mol, frac in mat.mpMolFrac.items():
                self.mpMolWt[mol] = self.mpMolWt.get(mol, 0) + frac * wt
        
        self.calcLoi()
        
    def calcLoi(self):
        wtFired = 0.0
        for wt in self.mpMolWt.values():
            wtFired += wt
        
        self.loi = self.wtTotal - wtFired
        
        for mol,wt in self.mpMolWt.items():
            self.mpMolWtFired[mol] = (wt / wtFired) * 100.0

##    def copy(self):
##        glaze = Glaze([])
##        glaze.mats = self.mats[:]
##        glaze.mpMatWt = self.mpMatWt.copy()
##        glaze.wtTotal = self.wtTotal
##        glaze.mpMolWt = self.mpMolWt.copy()
##        glaze.loi = self.loi
##        glaze.mpMolWtFired = self.mpMolWtFired.copy()
##        return glaze
##

    def oxides(self):
        "Returns dict of oxide=>wt"
        # BB better name for this?

        return self.mpMolWt.copy()

    def materials(self):
        mats = []
        for mat in self.mats:
            sMat = mat.name
            if self.mpMatFAdd[mat]:
                sMat = "+" + sMat
            mats.append(sMat)
        return mats

    def ingreds(self):
        ingreds = []
        for mat in self.mats:
            name = mat.name
            if self.mpMatFAdd[mat]:
                name = "+" + name
            ingreds.append((name, self.mpMatWt[mat]))
        return ingreds

    def grams(self, mol):
        return self.mpMolWt[mol]

    def moles(self, mol):
        return self.grams(mol) / mol.molWt

    @classmethod
    def fromBlend(cls, glazeA, fracA, glazeB, fracB, name):
        mpNameWt = {}
        for (sMat,wt) in glazeA.ingreds():
            mpNameWt[sMat] = wt * fracA
            
        for (sMat,wt) in glazeB.ingreds():
            wt *= fracB
            if sMat in mpNameWt:
                wt += mpNameWt[sMat]
            mpNameWt[sMat] = wt
        
        ingreds=[]
        for (sMat,wt) in sorted(mpNameWt.items(), key=lambda item : item[1], reverse=True):
            ingreds.append((sMat,wt))
            
        return Glaze(ingreds, name)

    @classmethod
    def fromGlazeChemFile(cls, filename):
        f = open(filename, "r")
        ingreds = []
        name = None
        for line in f.read().split('\n'):
            if len(line) == 0:
                continue
            
            match = re.match(r"""
                        (?P<field>[a-z_ ]+[a-z_])\s+=\s+
                        (?P<val>.*)""",
                        line, re.VERBOSE)
            if not match:
                print("can't parse", line)
                exit

            sField = match.group("field")
            sVal = match.group("val")

            if sField == "name":
                name = sVal
            elif sField == "component":
                sMat = sVal
            elif sField == "amount":
                ingreds.append((sMat, float(sVal)))
            elif sField == "addition":
                sMat = sVal
            elif sField == "addamount":
                ingreds.append(("+" + sMat, float(sVal)))

        return Glaze(ingreds, name)

    @classmethod
    def fromFile(cls, filename):
        if filename[-4:] == ".txt":
            return cls.fromGlazeChemFile(filename)
        else:
            print("unknown file type", filename[-4:])
            exit

    @classmethod          
    def fromFormula(cls, oxidesTarget, materials=None, limit=0.1, name=None):
        "Regenerate recipe from fired formula"

        mpMatFAdd = {}

        if materials is None:
            # by default, use full set in Materials
            # BB this duplicates things like CuO and CuCO3
            setMat = set()
            for s in Materials.values():
                s.add(s)
            materials = [mat.name for mat in setMat]
        
        oxides = {oxide:wt for oxide,wt in oxidesTarget.items() if wt>=limit}
        # Include oxides from unused materials--must penalize extra oxides!
        for sMat in materials:
            fAdd = (sMat[0] == '+')
            if fAdd:
                sMat = sMat[1:]
                
            mat = Materials[sMat]
            mpMatFAdd[mat] = fAdd
            
            for mol in mat.mpMolFrac.keys():
                if mol not in oxides:
                    oxides[mol] = 0.0

        # BB add extra row for combined fluxes (by molarity), scaled high,
        #  to ensure total flux equivalency?
        # BB score similarity by molarity instead of weight
        
        while True:
            a = np.zeros((len(oxides),len(materials)))
            b = np.zeros(len(oxides))
            for iOx,mol in enumerate(oxides):
                for iMat,sMat in enumerate(materials):
                    mat = Materials[sMat]
                    a[iOx,iMat] = mat.mpMolFrac.get(mol, 0)
                b[iOx] = oxides[mol]
            wts = np.linalg.lstsq(a,b,rcond=None)[0]
            if min(wts) > 0.0:
                break

            prunedMaterials = []
            for iMat,sMat in enumerate(materials):
                if wts[iMat] > 0.0:
                    prunedMaterials.append(sMat)
            materials = prunedMaterials

        ingreds = []
        for iMat,sMat in enumerate(materials):
            mat = Materials[sMat]
            sMat = mat.name
            if mpMatFAdd[mat]:
                sMat = "+" + sMat
            ingreds.append((sMat, wts[iMat]))

        ingreds.sort(key=lambda ingred : ingred[1], reverse=True)

        return cls(ingreds, name)

    def substitute(self, sMatOrig, *new):
        # BB always add +additives to end of list
        materials = self.materials()
        aMatNew = []
        for sMat in new:
            mat = Materials[sMat]
            aMatNew.append(mat.name)

        matOrig = Materials[sMatOrig].name
        i = materials.index(matOrig)
        materials = materials[0:i] + aMatNew + materials[i+1:]
        oxides = self.oxides()
        return Glaze.fromFormula(materials, oxides)

##    # BB How to do this? Forward optional arg list?
##    def sub(self, orig, *new):
##        return self.substitute(orig, new)

    def print(self, scale=1.0):
        print(self.name)

        for (sMat, wt) in self.ingreds():
            line = "{:<20}{:>7.1f}".format(sMat, wt*scale)
            print(line)
            
        line = "{:<20}{:>7.1f}".format("Total", self.wtTotal * scale)
        print(line)


    def printAnalysis(self):

        print(self.name)
        
        mols = []
        for mol in self.mpMolWt.keys():
            mols.append(mol)
        mols.sort()

        print("Percentage Analysis")
        line = " "*27
        for mol in mols:
            line += " {:>7}".format(mol.name)
        print(line + "{:>7}".format("LOI"))

        for mat in self.mats:
            wt = self.mpMatWt[mat]
            name = mat.name
            if self.mpMatFAdd[mat]:
                name = "+" + name
            line = "{:<20}{:>7.1f}".format(name, wt)
            for mol in mols:
                if mol in mat.mpMolFrac:
                    line += " {:>7.2f}".format(mat.mpMolFrac[mol] * wt)
                else:
                    line += " {:>7}".format("")
            print(line + "{:7.2f}".format(mat.loi * wt))

        line = "{:<20}{:>7.1f}".format("Total", self.wtTotal)
        for mol in mols:
            line += " {:>7.2f}".format(self.mpMolWt[mol])
        print(line + "{:7.2f}".format(self.loi))
        
        line = "{:<20}{:>7}".format("Fired %", "")
        for mol in mols:
            line += " {:>7.2f}".format(self.mpMolWtFired[mol])
        print(line)
        print()

        self.printSeger()

    def printSeger(self):
        
        # classify ingredients according to type
        fluxes = []
        glasses = []
        stabilizers = []

        for mol in self.mpMolWt:
            if mol.name == "SiO2" or mol.name == "TiO2":
                glasses.append(mol)
            elif mol.name == "Al2O3" or mol.name == "B2O3":
                stabilizers.append(mol)
            elif (mol.name == "PbO" or
                  mol.name == "Na2O" or
                  mol.name == "K2O" or
                  mol.name == "LiO" or
                  mol.name == "SrO" or
                  mol.name == "BaO" or
                  mol.name == "ZnO" or
                  mol.name == "CaO" or
                  mol.name == "MgO"):
                fluxes.append(mol)
            else:
                None # colorants not included in Seger formula

        sumFlux = 0.0
        for mol in fluxes:
            sumFlux += self.moles(mol)

        sumGlasses = 0.0
        for mol in glasses:
            sumGlasses += self.moles(mol)

        sumStab = 0.0
        for mol in stabilizers:
            sumStab += self.moles(mol)

        mpMolC = {}
        for mol,wt in self.mpMolWt.items():
            mpMolC[mol] = wt / (mol.molWt * sumFlux)

        print("{:8}{:>6}{:8}{:>6}{:8}{:>6}".format("", "RO", "", "R2O3", "", "RO2"))

        for i in range(max(len(fluxes), len(glasses), len(stabilizers))):
            line = ""
            if i < len(fluxes):
                line += "{:>8.2f}{:>6}".format(mpMolC[fluxes[i]], fluxes[i].name)
            else:
                line += "{:>8}{:>6}".format("", "")
            if i < len(stabilizers):
                line += "{:>8.2f}{:>6}".format(mpMolC[stabilizers[i]], stabilizers[i].name)
            else:
                line += "{:>8}{:>6}".format("", "")
            if i < len(glasses):
                line += "{:>8.2f}{:>6}".format(mpMolC[glasses[i]], glasses[i].name)
            else:
                line += "{:>8}{:>6}".format("", "")
            print(line)
        print()

        # R2O : RO 0.14 : 0.86

        # Silica:Alumina ratio predicts glossiness (higher is glossier)
        print("SiO2 : Al2O3 {:8.2f}".format(sumGlasses / sumStab))
        print()
    

