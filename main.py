# Add the import statements for functions from string_utils.py and equation_utils.py here
#import string_utils
import sympy
from sympy import symbols, Eq, solve as sympy_solve

def split_by_capitals(formula): #Splits a chemical formula (string) into a *list* of elements based on uppercase letters
  ret = []
  start = 0
  if not formula:  
    return []
  for i in range(1,len(formula)): 
    if formula[i].isupper():
      ret.append(formula[start:i])
      start = i
  ret.append(formula[start:])
  return ret

def split_at_digit(formula): #Splits an element string into a tuple of (element name, count)
  num = 0
  for i in range(len(formula)):
    if formula[i].isdigit():
      num=i
      break
  if num>0:
    return (formula[:i], int(formula[i:]))
  else:
    return (formula, 1)


def count_atoms_in_molecule(molecular_formula):
    """Takes a molecular formula (string) and returns a dictionary of atom counts.
    Example: 'H2O' → {'H': 2, 'O': 1}"""

    retDict = {} # Step 1: Initialize an empty dictionary to store atom counts

    for atom in split_by_capitals(molecular_formula):
        atom_name, atom_count = split_at_digit(atom)    # Step 2: Update the dictionary with the atom name and count
        retDict[atom_name]=atom_count

    return retDict # Step 3: Return the completed dictionary



def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count


#import equation_utils
# Add the import statements for necessary sympy functions here


ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
]

def generate_equation_for_element(compounds, coefficients, element):
    """Generates a symbolic equation for the given element from compounds and coefficients.  
    Example: For H in reactants [{'H': 2}, {'O': 4, 'H': 1}], coefficients [a0, a1], returns 2*a0 + a1."""
    equation = 0
    for i, compound in enumerate(compounds):
        if element in compound:
            equation += coefficients[i] * compound[element]
    return equation


def build_equations(reactant_atoms, product_atoms):
    """Builds a list of symbolic equations for each element to balance a chemical reaction.  
    Example: For H2 + O2 -> H2O, returns equations [2*a0 - 2*b0, a1 - b0]."""
    ## coefficients ##
    reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
    product_coefficients = list(symbols(f'b0:{len(product_atoms)}')) 
    product_coefficients = product_coefficients[:-1] + [1] # Ensure the last coefficient is 1

    ## equations ##
    equations = []
    for element in ELEMENTS:
        lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
        rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
        if lhs != 0 or rhs != 0:
            equations.append(Eq(lhs, rhs))

    return equations, reactant_coefficients + product_coefficients[:-1]


def my_solve(equations, coefficients):
    """Solves the system of equations for the coefficients of the reaction.  
    Example: For equations [2*a0 - 2*b0, a1 - b0], returns [1.0, 1.0]."""
    solution = sympy_solve(equations, coefficients)

    if len(solution) == len(coefficients):
        coefficient_values = list()
        for coefficient in coefficients:
            coefficient_values.append(float(solution[coefficient]))
        return coefficient_values

def balance_reaction(reaction): #"Fe2O3 + H2 -> Fe + H2O"

    # 1.parse reaction
    reactants, products = parse_chemical_reaction(reaction) # [""Fe2O3", "H2"], ["Fe", "H2O""]
    reactant_atoms = count_atoms_in_reaction(reactants) # [{"Fe":2, "O":1}, {"H":2}]
    product_atoms = count_atoms_in_reaction(products)

    # 2.build equation and solve
    equations, coefficients = build_equations(reactant_atoms, product_atoms)
    coefficients = my_solve(equations, coefficients) + [1]

    return coefficients # [1/3, 1, 2/3, 1]

