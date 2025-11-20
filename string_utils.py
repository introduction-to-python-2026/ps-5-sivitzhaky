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
