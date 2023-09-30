import numpy as np
from scipy.linalg import expm, norm
from math import radians

vdwr = {
    "H": 1.10,
    "O": 1.52,
    "N": 1.55,
    "C": 1.70,
    "S": 1.80
}

resname_3to1 = {
    "GLY": "G",
    "LEU": "L",
    "TYR": "Y",
    "SER": "S",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D",
    "ASN": "N",
    "PHE": "F",
    "ALA": "A",
    "LYS": "K",
    "ARG": "R",
    "HIS": "H",
    "CYS": "C",
    "VAL": "V",
    "PRO": "P",
    "TRP": "W",
    "ILE": "I",
    "MET": "M",
    "THR": "T",

    # Not standard
    "ASH": "D",
    "HSP": "H",
    "HSE": "H",
    "GLH": "E",
}


class Atom:
    def __init__(self):
        self.DataType = 'ATOM'  # "ATOM"/"HETATM"
        self.Name = ''  # Atom name
        self.AltLoc = ''  # Alternate location indicator.
        self.ResName = ''  # Residue name
        self.ChainId = 'A'  # Chain identifier
        self.ResSeq = 0  # Residue sequence number
        self.ResCode = ''  # Code for insertions of residues
        self.Coordin = np.array([0, 0, 0])  # (X,Y,Z) orthogonal
        self.Occup = 0.0  # Occupancy
        self.TempFac = 0.0  # Temperature factor
        self.Element = 'Xx'  # Element symbol
        self.Charge = 0.0  # Atom charge


class ZAtom:
    def __init__(self):
        self.Serial = 0
        self.Con = 0
        self.Dist = 0.0
        self.ValCon = 0
        self.ValAngle = 0.0
        self.DihCon = 0
        self.DihAngle = 0.0
        self.Name = 'Xx'


def M(axis, theta):
    return expm(np.cross(np.eye(3), axis / norm(axis) * theta))


def val_angle(vec1, vec2):
    dist01 = np.linalg.norm(vec1)
    dist02 = np.linalg.norm(vec2)
    dist12 = np.linalg.norm(vec1 - vec2)
    angle = (dist01 ** 2 + dist02 ** 2 - dist12 ** 2) / (2 * dist01 * dist02)
    return np.arccos(angle)


def read_pdb(fname):
    with open(fname, "r") as file:
        molecule = dict()
        for line in file.readlines():
            data_type = line[0:6].strip()
            if data_type not in ['ATOM', 'HETATM']:
                continue
            atom = Atom()
            atom.DataType = line[0:6].strip()
            num = int(line[6:11])
            atom.Name = line[12:16].strip()
            atom.AltLoc = line[16].strip()
            atom.ResName = line[17:20].strip()
            atom.ChainId = line[21].strip()
            atom.ResSeq = int(line[22:26])
            atom.ResCode = line[26].strip()
            atom.Coordin = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
            atom.Occup = 0.0  # float(line[54:60])
            atom.Tempfac = 0.0  # float(line[60:66])
            atom.Element = atom.Name[0]  # line[76:78].strip()

            molecule[num] = atom

    return molecule


def write_pdb(molecule, fname):
    with open(fname, "w") as f:
        for (idx, atom) in molecule.items():
            line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
                   '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}' \
                   '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
                a0=atom.DataType, a1=idx, a2=atom.Name, a3=atom.AltLoc,
                a4=atom.ResName, a5=atom.ChainId, a6=atom.ResSeq,
                a7=atom.ResCode, a8=atom.Coordin, a9=atom.Occup,
                a10=atom.TempFac, a11=atom.Element, s=' '
            )
            f.write(line)


def read_zmatrix(filename):
    zmatrix = list()
    with open(filename, "r") as f:
        resseq = int(f.readline().split()[-1])
        resname = f.readline().split()[-1]
        for line in f.readlines():
            words = line.split()
            zatom = ZAtom()
            zatom.Serial = int(words[0])
            zatom.Con = int(words[1])
            zatom.Dist = float(words[2])
            zatom.ValCon = int(words[3])
            zatom.ValAngle = radians(float(words[4]))
            zatom.DihCon = int(words[5])
            zatom.DihAngle = radians(float(words[6]))
            zatom.Name = words[7].strip()

            zmatrix.append(zatom)

    return resseq, resname, zmatrix


def read_sequence(fname):
    sequence = ""
    molecule = read_pdb(fname)
    for a in molecule.values():
        if a.Name == "CA":
            sequence += resname_3to1[a.ResName]
    return sequence


def read_coordinates(fname):
    coords = []
    molecule = read_pdb(fname)
    for atom in molecule.values():
        if atom.Name == "CA":
            coords.append(atom.Coordin)
    return coords
