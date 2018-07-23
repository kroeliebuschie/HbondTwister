#!/usr/bin/env python3
## 
## Author: Mehdi Nellen, Utrecht 2017

import sys
import math
import itertools
from copy import deepcopy
import numpy as np
import re


EXCLUDEDWATERS = [] # List containing waters that have to be excluded from the analysis
MINHBDIST = 3.5
MINHBANGLE = 2.1

class Atom:
    ''' class for atoms with element '''
    def __init__(self, atno, atom, resn, resi, xyzc, elem, conn, n, line):
        self.atno = atno # Unique atom number
        self.atom = atom # atom name 
        self.resn = resn # Residue/molecule name
        self.resi = resi # Residue/molecule integer/sequence number
        self.xyzc = xyzc # Coordinates xyz
        self.elem = elem # Element and charge
        self.conn = conn # List of connected elements
        self.hbty = ""   # participating in a hyrogen bond 
        self.hbon = 0    # hbond object
        self.n    = n    # line number in the pdb
        self.line = line # the raw pdb line (used for generating pdb)

    def calcDist(self, atom2):
        '''calculate distance between atoms''' 
        return np.linalg.norm(self.xyzc-atom2.xyzc)

    def calcAngle(self, atom1, atom3):
        '''calculates the angle located at atom2
        between the vectors atom1-atom2 and 
        atom2-atom3'''
        atom2 = self
        a21 = atom2.xyzc - atom1.xyzc
        a23 = atom2.xyzc - atom3.xyzc
        cos3 = np.dot(a21, a23) / (np.linalg.norm(a21) * np.linalg.norm(a23))
        angle = np.arccos(cos3)
        return angle
 
def findHbond(waters, molecules):
    '''This function searches for hydrogen bonds and 
    marks waters that form a hydrogen bond by putting a
    hydrogenbond object list in the water object'''
    # Waters as acceptors
    for water in waters:
        waterO = water.oxi
        #find Nitrogen Oxigen Fluoride
        for mol in molecules:
            mol.searchDonorHydrogens()
            ## Waters as acceptors
            for pair in mol.donorPairs:
                # search closest water oxigen
                if pair.checkHbond(waterO):
                    waterO.hbon = HydrogenBond(water, pair, "OA")
                    waterO.hbty = "OA"
                    water.oxi = waterO
                    water.hydrogenbonds += [HydrogenBond(pair, water, "OA")]

    #Waters as donors
    for mol in molecules:
        for acceptor in mol.Hacceptors:
            for water in waters:
                newhydros = set()
                for hydro in water.hydros:
                    HbDist = acceptor.calcDist(hydro)
                    if(MINHBDIST > HbDist):
                        if(MINHBANGLE < hydro.calcAngle(water.oxi, acceptor)):
                            hydro.hbon = HydrogenBond(acceptor, hydro, "HD")
                            hydro.hbty = "HD"
                            newhydros.add(hydro)
                            water.hydrogenbonds += [HydrogenBond(acceptor, hydro, "HD")]
                            print("#"*30)
                            print("HYDROGENBOND FOUND:")
                            print("angle:\t\t", 180*hydro.calcAngle(water.oxi, acceptor)/3.14)
                            print("dist:\t\t", HbDist)
                            print("acceptor:\t", acceptor.resi, " ", acceptor.atom)
                            print("donor:\t\t", hydro.resi, " ", hydro.atom)
                        else:
                            newhydros.add(hydro)
                    else:
                        newhydros.add(hydro)
                    water.hydros = newhydros

class HydrogenBond:
    def __init__(self, donor, acceptor, hbtype):
        ''' this class is instantiated with water molecule object,
        the hydrogen donor, hydrogen acceptor, and the type.
        There are 2 types of hydrogen bonds: 
        water oxigen as acceptor (OA) or water hydrogen as donor (HD)
        This is important for determining which hydrogens have to be rotated'''
        self.donor = donor
        self.acceptor = acceptor
        self.hbtype = hbtype

class DonorPair:
    '''class for a set of connected atoms which may form hydrogen bonds'''
    def __init__(self, hydrogen, NOF):
        self.hydrogen = hydrogen
        self.NOF = NOF

    def checkHbond(self, acceptorNOF):
        ''' this function checks whether there is an hydrogenbond formed between the 
        argument "acceptorNOF" and the donorpair. this is done by calculating the 
        angle and distance which have to match the restraints set in the global
        vairables'''
        # search closest water oxigen
        HbDist = self.hydrogen.calcDist(acceptorNOF) # calculate distance from acceptorNOF to donor H
        if(MINHBDIST > HbDist): # check for the minimum distance 
            if(MINHBANGLE < self.hydrogen.calcAngle(self.NOF, acceptorNOF)): # get the Hbond angle
                smallestDist = HbDist 
                print("#"*30)
                print("HYDROGENBOND FOUND:")
                print("angle:\t\t", 180*self.hydrogen.calcAngle(self.NOF, acceptorNOF)/3.14)
                print("dist:\t\t", HbDist)
                print("donor:\t\t", self.hydrogen.resi, self.hydrogen.atom)
                print("donorNOF:\t", self.NOF.resi, self.NOF.atom)
                print("acceptor:\t", acceptorNOF.resi, acceptorNOF.atom)
                return True
        # this return function replaces a bunch of 'else' statements
        return False

class Molecule:
    def __init__(self, num, *atoms):
        self.num = num           # pdb resi of the molecule
        self.atoms = set(atoms)
        self.Hacceptors = {}     # These are the NOFs
        self.donorHydrogens = [] # hydrogens connected to NOFs
        self.donorNOFs = []
        self.hydrogens = {}

    def addAtom(self, atom):
        '''adds atom to self.atoms'''
        self.atoms.add(atom)
    
    def filterElement(self, element):
        '''returns a set that matches the element which is being searched for'''
        return set(filter(lambda x: x.elem.startswith(element), self.atoms))

    def filterCon(self, c):
        '''returns a set that matches the element which is being searched for'''
        return set(filter(lambda x: x.atno == c, self.atoms))

    def filterNOF(self):
        '''returns a set that matches the elements NOF'''
        # Define regular expressions for Nitrogen Oxigen and Fluoride 
        # Sulfur can be added by simply adding S in the regex
        NOF = re.compile("[NOF][^A-Za-z]*")
        #filter for NOF elements
        return set(filter(lambda x: re.match(NOF, x.elem) is not None, self.atoms))

    def getHydrogens(self):
        return set(filter(lambda x: x.elem.startswith("H"), self.atoms))

    def searchAcceptors(self):
        self.Hacceptors = self.filterNOF()

    def searchDonorHydrogens(self):
        '''Search Hydrogen bond donors by looking for hydrogens 
           that are connected to a hydrogenbond acceptor '''
        if self.Hacceptors == {}:
            self.searchAcceptors()
        acceptorNo = [ a.atno for a in self.Hacceptors ]
        if self.hydrogens == {}:
            self.hydrogens = self.getHydrogens()
        self.donorPairs = [ DonorPair(h, self.filterCon(h.conn[0]).pop()) for h in self.hydrogens if h.conn[0] in acceptorNo ] # list of donor pairs
        self.donorHydrogens = [ h for h in self.hydrogens if h.conn[0] in acceptorNo ]
        
class Water(Molecule):
    ''' Water is a subclass of Molecule and conformers 
    can be created by twisting and tilting the water with
    respect to the bond'''
    def __init__(self, num, *atoms):
        Molecule.__init__(self, num, *atoms)
        self.oxi = self.filterElement("O").pop()
        self.hydros = self.filterElement("H")
        self.hydrogenbonds = []
        self.positions = []

    def rotate(self, rotateme, rotax, angle):
        '''Return the rotation matrix associated with counterclockwise rotation about
        the given axis by angle radians.'''
        rotax = np.asarray(rotax)
        rotax = rotax/math.sqrt(np.dot(rotax, rotax))
        a = math.cos(angle/2.0)
        b, c, d = -rotax*math.sin(angle/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rotmat = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                           [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                           [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        return np.dot(rotmat, rotateme )

    def findbinder(self):
        '''This function checks whether the atoms
        in the water molecule are marked as hydrogendonor 
        or hydrogen acceptor, if the water has 2 hydrogen 
        bonding atoms it has 0 degrees of freedom and None
        is returned'''
        count = 0
        for hydro in self.hydros:
            if hydro.hbty == "HD":
                count += 1
                binder = hydro
        if self.oxi.hbty == "OA":
            count += 1
            binder = self.oxi
        if count == 1:
            return binder
        else:
            return None

    def createWaterRotations(self, nrot):
        self.positions = [self]
        binder = self.findbinder()
        if binder is not None:
            print(binder.hbty)
            if binder.hbty == "HD":
                origin = self.oxi.xyzc
                rotax = binder.xyzc - origin
                for hydro in self.hydros:
                    if hydro.hbty == "":
                        rotateme = hydro.xyzc - origin
                        #self.positions = [ hydro ]
                        step = (2*np.pi)/nrot
                        for nstep in range(nrot)[1:]:
                            rotated = self.rotate(rotateme, rotax, step*nstep) #TODO step part can be left out if hydro is overwritten everytime
                            hydrocopy = deepcopy(hydro)
                            hydrocopy.xyzc = rotated + origin
                            self.positions += [Water(self.num, self.oxi, hydrocopy, binder)]
                            #[ print(hydro.xyzc) for water in self.positions for hydro in water.hydros]
            elif binder.hbty == "OA":
                origin = self.oxi.xyzc
                rotax =  self.hydrogenbonds[0].donor.hydrogen.xyzc - origin
                step = (2*np.pi)/nrot
                for nstep in range(nrot)[1:]:
                    hlist = [] # store rotaded hydrogens before they are added to the new water
                    for hydro in self.hydros:
                        rotateme = hydro.xyzc - origin
                        rotated = self.rotate(rotateme, rotax, step*nstep) #TODO step part can be left out if hydro is overwritten everytime
                        hydrocopy = deepcopy(hydro)
                        hydrocopy.xyzc = rotated + origin
                        hlist += [hydrocopy]
                    self.positions += [Water(self.num, self.oxi, hlist[0], hlist[1])]
            #[ print(hydro.xyzc) for water in self.positions for hydro in water.hydros]
        return

    def tilt():
        pass

class PDB:
    '''read file to a PDB object which is a collection of atoms'''
    def __init__(self, file_n):
        self.atoms = {} 
        self.remainingLINES = {}
        for n, line in enumerate(open(file_n)):
            if line.startswith("HETATM"):
                #self.atoms.add(
                self.atoms[int(line[6:11])] = Atom(int(line[6:11]),         # ATNO
                                              line[12:16].replace(" ", ""), # ATOM
                                              line[17:20].replace(" ", ""), # RESN
                                              int(line[22:26]),             # RESI 
                                              np.array((float(line[30:38]),  
                                                        float(line[38:46]),  
                                                        float(line[46:54]))),# XYZC
                                              line[77:80].replace(" ", ""), # ELEM
                                              [],                           # CONN
                                              n,                            # line number 
                                              line)                         # raw line 
                if line[17:20].replace(" ", "") != "HOH":
                    self.remainingLINES[n] = line
            elif line.startswith("CONECT"):
                self.remainingLINES[n] = line
                line = line.strip()[6:] # Remove trailing spaces and record name 
                # Add all connections of an atom
                self.atoms[int(line[0:6])].conn = [ int(line[i:i+6]) for i in range(6, len(line), 5) ]
            else:
                self.remainingLINES[n] = line
    def _toLine(self, atom):
        '''this function is able to convert molecules to 
        PDB text'''
        pdbline = atom.line
        pdbline = pdbline[:30] + "".join([ str(round(coord, 3)).rjust(8) for coord in atom.xyzc]) + pdbline[54:]
        return pdbline

    def writePDB(self, waters):
        ''' makes a pdb file out of the molecule and waters '''
        waterlists = [ water.positions for water in waters ]
        watercombos = list(itertools.product(*waterlists))
        for n, combo in enumerate(watercombos):
            dictcopy = deepcopy(self.remainingLINES)
            for water in combo:
                for atom in water.atoms:
                    dictcopy[atom.n] = self._toLine(atom)
            with open("testme%s.pdb" % str(n).zfill(4), "w") as output:
                l = 0 
                while len(dictcopy) > 0:
                    output.write(dictcopy.pop(l, None))
                    l += 1

    def filterAtoms(self, col, var, val, negate=False):
        ''' find atoms with value "val" for variable 
        "var" in set "col" and returns a filter object'''
        # Map dictionairy items (Atoms) and filter them based on var&val
        if negate:
            return filter(lambda x: getattr(x, var) != val, col)
        else:
            return filter(lambda x: getattr(x, var) == val, col)

    def getWaters(self):
        ''' create water object from atoms '''
        col = list(map(lambda a: a[1], self.atoms.items()))  # make a list out of all atoms
        allWaters = self.filterAtoms(col, "resn", "HOH")     # get water atoms by looking for HOH in resn
        waterNums = map(lambda w: w.resi, allWaters)         # get water residue numbers to group the atoms
        waterNums = list(set(waterNums))                     # filter for unique #TODO: use np.unique instead
        waters = { Water(num, *self.filterAtoms(col, "resi", num)) for num in waterNums }
        return waters

    def getMolecules(self):
        '''returns a dictionary with all non water molecules 
        as Molecule objects'''
        #TODO: Clone this method from waters
        col = list(map(lambda a: a[1], self.atoms.items()))    # make a list out of all atoms
        allMols = self.filterAtoms(col, "resn", "HOH", True)   # get molecules, the True argument inverts the selection
        molNums = map(lambda w: w.resi, allMols)               # get molecule residue numbers to group the atoms
        molNums = list(set(molNums))                           # filter for unique #TODO: use np.unique instead
        molecules = { Molecule(num, *self.filterAtoms(col, "resi", num)) for num in molNums }
        return molecules

if __name__ == "__main__":
    pdb = PDB(sys.argv[1])
    waters = pdb.getWaters()       # Dictionairy of water molecules 
    molecules = pdb.getMolecules() # Dictionairy of non water molecules
    findHbond(waters, molecules)
    [ water.createWaterRotations(4) for water in waters ]
    pdb.writePDB(waters)
