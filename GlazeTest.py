from Glaze import *

##haynes.printAnalysis()
##haynes.printSeger()
##haynes2 = haynes.substitute("Nepheline Syenite", "Custer Feldspar", "Kona F-4", "+RIO")
##haynes2.printAnalysis()

##g = Glaze.fromFile("glazes/Glazy_ID_2462_GlazeChem.txt")
##g.printAnalysis()


cone_6_mat = Glaze([\
    ("Custer Feldspar", 63.6),
    ("Whiting", 18.3),\
    ("Kaolin", 9.1),\
    ("Talc", 4.5),\
    ("Zinc Oxide", 4.5)\
    ])
#cone_6_mat.printAnalysis()

haynes = Glaze([\
    ("Nepheline Syenite", 45),
    ("Silica", 30),
    ("Dolomite", 10),
    ("Whiting", 8),
    ("Talc", 7),
    ("+Bentonite", 3)], "Haynes White")
#haynes.printAnalysis()

#haynes_olive = Glaze(haynes.ingreds() + [("+RIO", 2)], "Haynes Olive")
#haynes_olive.printAnalysis()

##materials = haynes.materials()
##materials.append("EPK")
#materials = ["EPK", "Custer Feldspar", "+Bentonite"]
#oxides = haynes.oxides()                           
#replaceOxide(oxides, "Na2O", 0.5, "MgO")
#mulOxide(oxides, "SiO2", 1.1)
#new_haynes = Glaze.recalc(materials, oxides) 
#new_haynes.printAnalysis()

##test = Glaze([("EPK", 100)])
##test.printAnalysis()
##test2 = Glaze.fromFormula(["Nepheline Syenite", "Custer Feldspar", "Barnard Clay", "Talc", "Silica", "Whiting", "+Bentonite", "+RIO", "+Rutile"], test.oxides())
##test2.printAnalysis()


bto = Glaze([\
    ("Custer Feldspar", 490),
    ("Cornish Stone", 180),
    ("Whiting", 180),
    ("EPK", 60),
    ("OM-4", 90),
    ("Zinc Oxide", 30),
    ("+RIO", 60),
    ("+Rutile", 40)
    ], "B Temple Orange")
#bto.printAnalysis()

#bto2 = bto.substitute("Cornish Stone", "Gerstley Borate")
#bto2.printAnalysis() # virtually no GB

#bto3 = bto.substitute("Cornish Stone")
#bto3.printAnalysis() # removes EPK entirely! Too much (20%) ball clay?

#bto4 = bto.substitute("Cornish Stone", "Silica")
#bto4.printAnalysis()
#bto4.print(5)


willie_helix_2 = Glaze([
    ("Nepheline Syenite", 40),
    ("Whiting", 19),
    ("Silica", 30),
    ("EPK", 11),
    ("+CuO", 5)
    ], "Willie Helix 2")
#willie_helix_2.printAnalysis()


yellow_salt = Glaze([
    ("Nepheline Syenite", 71.6),
    ("OM-4", 4.8),
    ("Dolomite", 23.6),
    ("Zircopax", 17.9),
    ("+Bentonite", 4.0),
    ("+RIO", 1.1)
    ], "Yellow Salt")
#yellow_salt.printAnalysis()


sage1=Glaze.fromBlend(yellow_salt, .80, willie_helix_2, .20, "ys80wh20")
#sage1.printAnalysis()

#sage2=Glaze.fromFormula(sage1.oxides())
#sage2.printAnalysis()

sage3=Glaze.fromFormula(sage1.oxides(), sage1.materials(), name="Sage3(reform)")
#sage3.printAnalysis()

sage4=Glaze(sage3.ingreds() + [("+Bentonite", 2.0)], name="Sage")
#sage4.printAnalysis()
#sage4.print(scale = 5)

#interesting to try: allow MgO/CaO ratio to vary, use only Dolomite

ashten=Glaze([
    ("Soda Spar", 5940),
    ("EPK", 360),
    ("Whiting", 720),
    ("Silica", 1080),
    ("Wood Ash", 900),
    ("+Yellow Ochre", 945)],
     "Wood Ash Tenmoku")

ashten.print(1)
ashten.print(5)

